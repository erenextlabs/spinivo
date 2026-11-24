"""Transform, reference and integrate 1D proton spectra on the acquired intensity scale."""

from dataclasses import dataclass, replace
from math import inf, isfinite, sqrt
from numbers import Real
from pathlib import Path

import numpy as np
from nmrglue.process.proc_autophase import autops
from nmrglue.process.proc_base import fft_positive
from numpy.typing import NDArray

from .filters import DelayDecision, remove_delay, resolve_delay
from .integral_checks import (
    BaselineRecovery, LockedBaselineRecovery, SpikeAnchorFlips, _BASELINE_SENSITIVITY_LIMIT, _measure_baseline_recovery,
)
from .readers.model import RawFID


@dataclass(frozen=True)
class ProcessingSettings:
    """Phase is p0 at the pivot, increasing by p1 over one sweep left to right.

    A missing pivot uses the referenced carrier. ACME uses p0 and p1 as initial
    values and sets positive polarity at the largest magnitude point. Lineshape
    assumes positive peaks and estimates both phase values from the spectrum.
    Set phase_lb_hz to estimate phase with separate broadening from lb_hz.
    Set baseline_bridge to join nearby anchors across declared regions on the smooth branch.
    Set lock_spike_anchors to hold the nominal anchors and branch fixed in spike trials.
    """

    lb_hz: float = 0.0
    zero_fill_factor: int = 2
    reference_shift_ppm: float = 0.0
    reference_source: str | None = None
    p0_deg: float = 0.0
    p1_deg: float = 0.0
    pivot_ppm: float | None = None
    autophase_method: str | None = None
    baseline_windows_ppm: tuple[tuple[float, float], ...] = ()
    smooth_baseline: bool = False
    phase_lb_hz: float | None = None
    baseline_bridge: bool = False
    lock_spike_anchors: bool = False


# Estimate phase on a 1 Hz copy for the lineshape fit.
# Integrate without broadening to avoid line-width-dependent losses in fixed windows.
QUANTITATIVE_SETTINGS = ProcessingSettings(
    zero_fill_factor=2, lb_hz=0.0, phase_lb_hz=1.0, autophase_method="lineshape", smooth_baseline=True,
    baseline_bridge=True, lock_spike_anchors=True,
)


def uses_quantitative_settings(settings: ProcessingSettings) -> bool:
    """A documented reference shift is allowed within the quantitative profile."""
    source = settings.reference_source
    if settings.reference_shift_ppm and (not isinstance(source, str) or not source.strip()):
        return False
    windows = tuple(tuple(window) for window in settings.baseline_windows_ppm)
    return replace(settings, reference_shift_ppm=0.0, reference_source=None,
                   baseline_windows_ppm=windows) == QUANTITATIVE_SETTINGS


@dataclass(frozen=True)
class LineshapePeak:
    ppm: float
    half_width_ppm: float
    height: float
    phase_deg: float
    relative_rms: float


@dataclass(frozen=True)
class LineshapePhaseFit:
    """Phase estimation method, accepted peaks and fit quality."""

    peaks: tuple[LineshapePeak, ...]
    slope_deg: float
    phase_rms_deg: float | None
    method: str
    note: str | None = None


@dataclass(frozen=True)
class BaselineShapeAreas:
    """Signed areas in intensity*Hz under the three fixed baseline candidates."""

    production: float
    anchor_median: float
    bracketing_line: float | None


@dataclass(frozen=True)
class BaselineRegion:
    """Distances to the bracketing anchors in low to high ppm order."""

    bounds_ppm: tuple[float, float]
    anchor_distances_ppm: tuple[float | None, float | None]
    half_block_area: float | None
    half_block_relative_change: float | None
    recovery: BaselineRecovery | None = None
    bridged: bool = False
    shape_areas: BaselineShapeAreas | None = None


@dataclass(frozen=True)
class BaselineSelection:
    branch: str
    eligible_blocks: int
    adjacent_pairs: int
    median_residual: float | None
    residual_mad: float | None
    noise: float | None
    reason: str | None = None


@dataclass(frozen=True)
class SmoothBaselineFit:
    """Automatic baseline diagnostics. Constant-branch levels have fitted wings removed."""

    anchors: tuple[tuple[float, float], ...]  # Block center in ppm and intensity used for fitting.
    regions: tuple[BaselineRegion, ...]
    selection: BaselineSelection | None = None
    block_size_points: int | None = None
    half_block_selection: BaselineSelection | None = None


@dataclass(frozen=True)
class Spectrum:
    """Complex spectrum and effective settings. Baseline coefficients are intercept and slope per ppm."""

    data: NDArray[np.complex128]
    hz: NDArray[np.float64]
    ppm: NDArray[np.float64]
    source: Path
    intensity_note: str
    delay: DelayDecision
    delay_removal_method: str
    tail_trim_points: int
    acquired_points: int
    transformed_points: int
    dwell_s: float
    spectral_width_hz: float
    observe_frequency_mhz: float
    reference_ppm: float
    ppm_limits: tuple[float, float]
    settings: ProcessingSettings
    baseline_coefficients: tuple[float, float] | None
    first_point_weight: float
    transform: str
    lineshape_phase_fit: LineshapePhaseFit | None = None
    smooth_baseline_fit: SmoothBaselineFit | None = None
    linear_baseline_recoveries: tuple[BaselineRecovery, ...] = ()
    smooth_baseline_curve: NDArray[np.float64] | float | None = None  # Subtracted real baseline.
    baseline_limit: float = _BASELINE_SENSITIVITY_LIMIT


@dataclass(frozen=True)
class Integral:
    bounds_ppm: tuple[float, float]
    bounds_hz: tuple[float, float]
    area: float
    unit: str = "intensity*Hz"


def _ordered_bounds(bounds_ppm, ppm_limits):
    low, high = sorted(bounds_ppm)
    if not ppm_limits[0] <= low < high <= ppm_limits[1]:
        raise ValueError("ppm bounds must span a nonzero interval within the spectral width")
    return low, high


def _fit_baseline(data, ppm, windows, ppm_limits):
    anchors = np.zeros(ppm.shape, dtype=bool)
    for window in windows:
        low, high = _ordered_bounds(window, ppm_limits)
        in_window = (low <= ppm) & (ppm <= high)
        if not in_window.any():
            raise ValueError(f"Baseline window {window} contains no spectrum points")
        anchors |= in_window
    if np.count_nonzero(anchors) < 2:
        raise ValueError("Baseline windows must contain at least two distinct spectrum points")
    intercept, slope = np.polynomial.polynomial.polyfit(ppm[anchors], data.real[anchors], 1)
    data.real -= intercept + slope * ppm
    return float(intercept), float(slope)


def _fit_peak(ppm, data, index, half_width):
    # Fit over three apparent half-widths with at least five samples on each side.
    spacing = ppm[0] - ppm[1]
    radius = max(5, int(np.ceil(3 * half_width / spacing)))
    if index < radius or index + radius >= ppm.size:
        return None
    window = slice(index - radius, index + radius + 1)
    x = (ppm[window] - ppm[index]) / half_width
    signal = data[window] / abs(data[index])

    def fit(center, width):
        peak = 1 / (1 + 1j * (x - center) / width)
        design = np.column_stack((peak, np.ones_like(x), x))
        coefficients = np.linalg.lstsq(design, signal, rcond=None)[0]
        return np.linalg.norm(signal - design @ coefficients), coefficients[0], np.linalg.norm(peak)

    # Solve the complex linear terms exactly at each trial center and width.
    center, width, step = 0.0, 1.0, 0.25
    error, amplitude, norm = fit(center, width)
    for _ in range(32):
        start_center, start_width = center, width
        for dc, dw in ((-step, 0), (step, 0), (0, -step), (0, step)):
            trial_center, trial_width = start_center + dc, start_width + dw
            if abs(trial_center) > 1 or not 0.25 <= trial_width <= 2:
                continue
            trial_error, trial_amplitude, trial_norm = fit(trial_center, trial_width)
            if trial_error < error:
                center, width = trial_center, trial_width
                error, amplitude, norm = trial_error, trial_amplitude, trial_norm
        if (center, width) == (start_center, start_width):
            step /= 2
            if step < 0.0001:
                break
    height = abs(amplitude) * abs(data[index])
    if height == 0 or width * half_width < spacing / 2 or error > 0.05 * abs(amplitude) * norm:
        return None
    return LineshapePeak(
        float(ppm[index] + center * half_width), float(width * half_width), float(height),
        float(-np.angle(amplitude, deg=True)), float(error / (abs(amplitude) * norm)),
    )


def _wrapped(degrees):
    return (degrees + 180) % 360 - 180


def _peak_prominence(magnitude, index):
    height = magnitude[index]
    left = np.flatnonzero(magnitude[:index] > height)
    right = np.flatnonzero(magnitude[index + 1:] > height)
    start = left[-1] + 1 if left.size else 0
    stop = index + 1 + right[0] if right.size else magnitude.size
    return height - max(magnitude[start:index + 1].min(), magnitude[index:stop].min())


def _lineshape_peaks(data, ppm):
    if data.size < 3:
        return []
    magnitude = np.abs(data)
    threshold = max(0.002 * magnitude.max(), 12 * np.median(abs(np.diff(data))))
    candidates = np.flatnonzero(
        (magnitude[1:-1] >= magnitude[:-2]) & (magnitude[1:-1] > magnitude[2:])
    ) + 1
    prominences = {}

    def prominence(index):
        if index not in prominences:
            prominences[index] = _peak_prominence(magnitude, index)
        return prominences[index]

    spacing = ppm[0] - ppm[1]
    limit = max(5, ppm.size // 100)
    peaks = []
    for index in candidates[magnitude[candidates] > threshold]:
        if prominence(index) <= threshold:
            continue
        level = magnitude[index] / np.sqrt(2)
        left, right = index, index
        while left > max(0, index - limit) and magnitude[left] > level:
            left -= 1
        while right < min(ppm.size - 1, index + limit) and magnitude[right] > level:
            right += 1
        if magnitude[left] > level or magnitude[right] > level:
            continue
        half_width = (right - left) * spacing / 2
        radius = 1.5 * (right - left)
        start = candidates.searchsorted(index - radius)
        stop = candidates.searchsorted(index + radius, side="right")
        if any(other != index and magnitude[other] >= 0.1 * magnitude[index]
               and prominence(other) >= 0.1 * magnitude[index] for other in candidates[start:stop]):
            continue
        peak = _fit_peak(ppm, data, index, half_width)
        if peak is not None:
            peaks.append(peak)
    return peaks


def _lineshape_phase(data, ppm, pivot, width_ppm, edges):
    peaks = _lineshape_peaks(data, ppm)
    if not peaks:
        raise ValueError("Lineshape autophase found no accepted peaks")
    position = (pivot - np.array([peak.ppm for peak in peaks])) / width_ppm
    edge_position = (pivot - np.asarray(edges)) / width_ppm
    phases = np.array([peak.phase_deg for peak in peaks])
    height = np.array([peak.height for peak in peaks])
    error = np.array([peak.relative_rms for peak in peaks])
    weights = np.sqrt(height / height.max()) / np.maximum(error, 0.01) ** 2
    slopes = np.arange(-360, 361)
    detrended = phases - slopes[:, None] * position
    intercepts = np.angle(np.exp(1j * np.deg2rad(detrended)) @ weights, deg=True)
    errors = _wrapped(detrended - intercepts[:, None])
    scores = np.average(np.minimum(errors ** 2, 20 ** 2), axis=1, weights=weights)
    best = np.argmin(scores)
    inliers = abs(errors[best]) <= 20

    if np.count_nonzero(inliers) < 2 or np.ptp(position[inliers]) < 0.05:
        strongest = np.argmax(weights)
        if 360 * max(abs(edge_position - position[strongest])) <= 5:
            fit = LineshapePhaseFit((peaks[strongest],), 0.0, None, "single_peak", "slope not determined")
            return phases[strongest], 0.0, fit
        raise ValueError("Lineshape autophase has insufficient peak span for the requested regions")

    minima = (scores <= np.r_[np.inf, scores[:-1]]) & (scores <= np.r_[scores[1:], np.inf])
    rivals = minima & (abs(slopes - slopes[best]) > 30) & (scores <= scores[best] + 4)
    edge_difference = ((intercepts[rivals] - intercepts[best])[:, None]
                       + (slopes[rivals] - slopes[best])[:, None] * edge_position)
    if np.any(abs(_wrapped(edge_difference)) > 5):
        raise ValueError("Lineshape autophase has ambiguous phase at a requested region edge")

    design = np.column_stack((np.ones(np.count_nonzero(inliers)), position[inliers]))
    unwrapped = errors[best, inliers] + intercepts[best] + slopes[best] * position[inliers]
    root = np.sqrt(weights[inliers])
    ramp = np.linalg.lstsq(design * root[:, None], unwrapped * root, rcond=None)[0]
    rms = np.sqrt(np.average((unwrapped - design @ ramp) ** 2, weights=weights[inliers]))
    if rms > 10 or abs(ramp[1]) > 360:
        raise ValueError("Lineshape autophase exceeds the phase RMS or slope limit")
    accepted = tuple(peak for peak, keep in zip(peaks, inliers) if keep)
    return _wrapped(ramp[0]), ramp[1], LineshapePhaseFit(accepted, float(ramp[1]), float(rms), "ramp")


def _smooth_levels(levels, anchors, penalty):
    # Cholesky factorization of diag(anchors) + penalty * D2.T @ D2, in three bands.
    diagonal = anchors.astype(float)
    diagonal[:-2] += penalty
    diagonal[1:-1] += 4 * penalty
    diagonal[2:] += penalty
    first = np.zeros_like(diagonal)
    first[1:-1] -= 2 * penalty
    first[2:] -= 2 * penalty
    second = np.full_like(diagonal, penalty)
    second[:2] = 0
    result = np.where(anchors, levels, 0.0)
    for i in range(diagonal.size):
        if i >= 2:
            second[i] /= diagonal[i - 2]
            result[i] -= second[i] * result[i - 2]
        if i >= 1:
            first[i] = (first[i] - second[i] * first[i - 1]) / diagonal[i - 1]
            result[i] -= first[i] * result[i - 1]
        diagonal[i] = np.sqrt(diagonal[i] - first[i] ** 2 - second[i] ** 2)
        result[i] /= diagonal[i]
    for i in reversed(range(diagonal.size)):
        if i + 2 < diagonal.size:
            result[i] -= second[i + 2] * result[i + 2]
        if i + 1 < diagonal.size:
            result[i] -= first[i + 1] * result[i + 1]
        result[i] /= diagonal[i]
    return result


def _mad(values):
    return float(np.median(abs(values - np.median(values))))


def _select_baseline(centers, levels, anchors, phase_fit):
    eligible = anchors.copy()
    residuals = levels.copy()
    if phase_fit is not None:
        for peak in phase_fit.peaks:
            distance = (centers - peak.ppm) / peak.half_width_ppm
            eligible &= abs(distance) > 80  # 40 fitted full widths.
            residuals -= peak.height / (1 + distance ** 2)
    kept = residuals[eligible]
    adjacent = np.diff(residuals)[eligible[:-1] & eligible[1:]]
    median = float(np.median(kept)) if kept.size else None
    spread = _mad(kept) if kept.size else None
    noise = _mad(adjacent) / sqrt(2) if adjacent.size else None
    reason = None
    if adjacent.size < 11:
        reason = "Fewer than 11 adjacent eligible block pairs"
    elif spread > 2 * noise:
        reason = "Residual MAD exceeds twice the adjacent block noise"
    branch = "constant" if reason is None else "smooth"
    return BaselineSelection(branch, kept.size, adjacent.size, median, spread, noise, reason), eligible, residuals


@dataclass(frozen=True)
class _BaselineCurve:
    curve: NDArray[np.float64] | float
    anchor_ppm: NDArray[np.float64]
    levels: NDArray[np.float64]
    selection: BaselineSelection
    bridged: list[bool]
    placement_curve: NDArray[np.float64] | None
    anchor_mask: NDArray[np.bool_]
    added_ppm: tuple[float, ...]
    removed_ppm: tuple[float, ...]


def _baseline_curve(spectrum, regions, block_size, offset=0, *, measure_placement=False, anchor_lock=None):
    stop = offset + ((spectrum.ppm.size - offset) // block_size) * block_size
    blocks = spectrum.data.real[offset:stop].reshape(-1, block_size)
    if len(blocks) < 2:
        return None
    high = spectrum.ppm[offset:stop:block_size]
    low = spectrum.ppm[offset + block_size - 1:stop:block_size]
    centers = (high + low) / 2
    x = np.linspace(-0.5, 0.5, block_size)
    slope = blocks @ x / (x @ x)
    residual = blocks - blocks.mean(axis=1)[:, None] - slope[:, None] * x
    error = 1.4826 * np.median(abs(residual), axis=1)
    floor = np.percentile(error, 25)
    anchors = (error <= 3 * floor) & (abs(slope) <= max(6 * floor, 0.0001 * abs(spectrum.data).max()))
    for start, end in regions:
        anchors &= (high < start) | (low > end)
    added_ppm, removed_ppm = (), ()
    if anchor_lock is not None:
        mask, _ = anchor_lock
        # The trial's free mask is evidence, not an input to the locked fit.
        added_ppm = tuple(map(float, centers[anchors & ~mask]))
        removed_ppm = tuple(map(float, centers[mask & ~anchors]))
        anchors = mask
    if np.count_nonzero(anchors) < 2:
        return None
    levels = np.median(blocks, axis=1)
    # Fitted wings describe the quantitative trace only for a direct phase fit.
    phase_fit = spectrum.lineshape_phase_fit if spectrum.settings.phase_lb_hz is None else None
    selection, eligible, residuals = _select_baseline(centers, levels, anchors, phase_fit)
    if anchor_lock is not None and selection.branch != anchor_lock[1]:
        selection = replace(selection, branch=anchor_lock[1], reason="Nominal branch held fixed in spike trial")
    bridged = [False] * len(regions)
    if selection.branch == "constant":
        return _BaselineCurve(selection.median_residual, centers[eligible], residuals[eligible], selection, bridged,
                              None, anchors, added_ppm, removed_ppm)
    smooth = _smooth_levels(levels, anchors, (0.025 / (centers[0] - centers[1])) ** 4)
    straight = smooth.copy() if measure_placement else None
    if spectrum.settings.baseline_bridge or straight is not None:
        for index, (lo, hi) in enumerate(regions):
            below = np.flatnonzero(anchors & (centers < lo))
            above = np.flatnonzero(anchors & (centers > hi))
            if above.size and below.size:
                # The block centers descend in ppm.
                i, j = above[-1], below[0]
                close = max(centers[i] - hi, lo - centers[j]) <= (hi - lo) / 2
                use_bridge = spectrum.settings.baseline_bridge and close
                if not use_bridge and straight is None:
                    continue
                line = np.linspace(smooth[i], smooth[j], j - i + 1)
                if straight is not None:
                    straight[i:j + 1] = line
                if use_bridge:
                    smooth[i:j + 1] = line
                    bridged[index] = True
    curve = np.interp(spectrum.ppm, centers[::-1], smooth[::-1])
    placement_curve = np.interp(spectrum.ppm, centers[::-1], straight[::-1]) if straight is not None else None
    return _BaselineCurve(curve, centers[anchors], levels[anchors], selection, bridged, placement_curve,
                          anchors, added_ppm, removed_ppm)


def _baseline_shape_areas(spectrum, bounds, nominal, median_spectrum, fitted):
    """Evaluate the alternatives on the same phased trace and integration grid."""
    ppm, levels = fitted.anchor_ppm, fitted.levels
    below = np.flatnonzero(ppm < bounds[0])
    above = np.flatnonzero(ppm > bounds[1])
    line_area = None
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        median_area = integrate(median_spectrum, bounds).area
        if below.size and above.size:
            i, j = below[0], above[-1]  # Anchor centers descend in ppm.
            slope = (levels[j] - levels[i]) / (ppm[j] - ppm[i])
            line = levels[i] + slope * (spectrum.ppm - ppm[i])
            line_area = integrate(replace(spectrum, data=spectrum.data - line), bounds).area
    return BaselineShapeAreas(nominal, median_area, line_area)


def _apply_automatic_baseline(spectrum, regions, *, baseline_shape=False):
    block_size = max(8, round(4 / (spectrum.hz[0] - spectrum.hz[1]))) if spectrum.hz.size > 1 else 8
    fitted = _baseline_curve(spectrum, regions, block_size, measure_placement=bool(regions))
    if fitted is None:
        raise ValueError("Smooth baseline requires at least two quiet anchor blocks")
    curve, anchor_ppm, placement_curve = fitted.curve, fitted.anchor_ppm, fitted.placement_curve
    shifted = _baseline_curve(spectrum, regions, block_size, block_size // 2) if regions else None
    shifted_spectrum = replace(spectrum, data=spectrum.data - shifted.curve) if shifted is not None else None
    corrected = replace(spectrum, data=spectrum.data - curve)
    straight_spectrum = replace(spectrum, data=spectrum.data - placement_curve) if placement_curve is not None else None
    anchor_lock = (fitted.anchor_mask, fitted.selection.branch) if spectrum.settings.lock_spike_anchors else None
    median_spectrum = None
    if baseline_shape and regions:
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            median_spectrum = replace(spectrum, data=spectrum.data - np.median(fitted.levels))

    def refit(trial, position_percent):
        fitted_trial = _baseline_curve(trial, regions, block_size, anchor_lock=anchor_lock)
        if fitted_trial is None:
            raise ValueError("Added signal leaves fewer than two quiet anchor blocks")
        if fitted_trial.added_ppm or fitted_trial.removed_ppm:
            flips.append(SpikeAnchorFlips(position_percent, fitted_trial.added_ppm, fitted_trial.removed_ppm))
        return replace(trial, data=trial.data - fitted_trial.curve)

    diagnostics = []
    for bounds, used_bridge in zip(regions, fitted.bridged):
        low, high = bounds
        distances = (low - anchor_ppm[anchor_ppm < low].max(initial=-np.inf),
                     anchor_ppm[anchor_ppm > high].min(initial=np.inf) - high)
        nominal = integrate(corrected, bounds).area
        alternate = integrate(shifted_spectrum, bounds).area if shifted_spectrum is not None else None
        relative = (alternate - nominal) / nominal if alternate is not None and nominal else None
        straight_area = nominal if fitted.selection.branch == "constant" else None
        if straight_spectrum is not None and all(map(isfinite, distances)):
            straight_area = integrate(straight_spectrum, bounds).area
        flips = []
        recovery = _measure_baseline_recovery(spectrum, bounds, nominal, refit, straight_area,
                                              limit=spectrum.baseline_limit)
        if anchor_lock is not None:
            recovery = LockedBaselineRecovery(**vars(recovery), anchor_flips=tuple(flips))
        shape_areas = (_baseline_shape_areas(spectrum, bounds, nominal, median_spectrum, fitted)
                       if median_spectrum is not None else None)
        diagnostics.append(BaselineRegion(bounds, tuple(float(d) if isfinite(d) else None for d in distances),
                                          alternate, relative, recovery, used_bridge, shape_areas))
    anchors = tuple(zip(map(float, anchor_ppm), map(float, fitted.levels)))
    fit = SmoothBaselineFit(anchors, tuple(diagnostics), fitted.selection, block_size,
                            shifted.selection if shifted is not None else None)
    return replace(corrected, smooth_baseline_fit=fit, smooth_baseline_curve=curve)


def _apply_linear_baseline(spectrum, regions):
    def refit(trial, _position_percent=None):
        coefficients = _fit_baseline(trial.data, trial.ppm, trial.settings.baseline_windows_ppm, trial.ppm_limits)
        return replace(trial, baseline_coefficients=coefficients)

    corrected = refit(replace(spectrum, data=spectrum.data.copy()))
    recoveries = []
    for bounds in regions:
        nominal = integrate(corrected, bounds).area
        recoveries.append(_measure_baseline_recovery(spectrum, bounds, nominal, refit, nominal,
                                                     limit=spectrum.baseline_limit))
    return replace(corrected, linear_baseline_recoveries=tuple(recoveries))


def _transform_fid(fid, lb_hz, dwell, size, first_point_weight):
    """Apply the window and first-point weight in place before transforming."""
    if lb_hz:
        fid *= np.exp(-np.pi * lb_hz * dwell * np.arange(fid.size))
    fid[0] *= first_point_weight
    return fft_positive(np.pad(fid, (0, size - fid.size))) * dwell


def process(raw: RawFID, settings: ProcessingSettings = ProcessingSettings(), *,
            regions_ppm: tuple[tuple[float, float], ...] = (), baseline_limit: float = _BASELINE_SENSITIVITY_LIMIT) -> Spectrum:
    """Remove group delay, transform the FID and apply the requested corrections."""
    correction = remove_delay(raw, resolve_delay(raw))
    sw, observe, carrier = raw.spectral_width_hz, raw.observe_frequency_mhz, raw.carrier_ppm
    for name, value in (("spectral_width_hz", sw), ("observe_frequency_mhz", observe), ("carrier_ppm", carrier)):
        if value is None or not isfinite(value) or (name != "carrier_ppm" and value <= 0):
            raise ValueError(f"Missing or invalid {name}")
    if not isfinite(settings.lb_hz) or settings.lb_hz < 0:
        raise ValueError("lb_hz must be finite and nonnegative")
    if type(settings.zero_fill_factor) is not int or settings.zero_fill_factor < 1:
        raise ValueError("zero_fill_factor must be a positive integer")
    if not all(map(isfinite, (settings.reference_shift_ppm, settings.p0_deg, settings.p1_deg))):
        raise ValueError("Reference shift and phase values must be finite")
    if not isinstance(baseline_limit, Real) or isinstance(baseline_limit, bool) or not 0 <= baseline_limit < inf:
        raise ValueError("baseline_limit must be finite and nonnegative")
    if settings.autophase_method not in (None, "acme", "lineshape"):
        raise ValueError("autophase_method must be 'acme', 'lineshape' or None")
    phase_lb = settings.phase_lb_hz
    if phase_lb is not None:
        if not isfinite(phase_lb) or phase_lb < 0:
            raise ValueError("phase_lb_hz must be finite and nonnegative")
        if settings.autophase_method is None:
            raise ValueError("phase_lb_hz requires an autophase_method")
    if settings.smooth_baseline and settings.baseline_windows_ppm:
        raise ValueError("Choose either a smooth baseline or explicit linear baseline windows")
    if raw.complex_points < 1 or raw.fid.shape != (raw.complex_points,):
        raise ValueError("FID shape does not match the acquired complex points")
    if not np.isfinite(correction.fid).all():
        raise ValueError("Delay-corrected FID must contain only finite samples")

    dwell = 1 / sw
    size = raw.complex_points * settings.zero_fill_factor
    first_point_weight = 0.5
    phase_data = None
    if phase_lb is not None:
        phase_data = _transform_fid(correction.fid.copy(), phase_lb, dwell, size, first_point_weight)
    data = _transform_fid(correction.fid, settings.lb_hz, dwell, size, first_point_weight)
    if phase_data is None:
        phase_data = data

    # Both vendor fixtures match their processed peaks when fft_positive bins
    # are assigned the opposite frequencies. fftshift then gives descending Hz.
    hz = (size // 2 - np.arange(size)) / size * sw
    reference = carrier + settings.reference_shift_ppm
    width_ppm = sw / observe
    ppm = reference + hz / observe
    ppm_limits = (reference - width_ppm / 2, reference + width_ppm / 2)
    pivot = reference if settings.pivot_ppm is None else settings.pivot_ppm
    if not ppm_limits[0] <= pivot <= ppm_limits[1]:
        raise ValueError("pivot_ppm must lie within the referenced spectral width")
    regions = tuple(_ordered_bounds(bounds, (ppm[-1], ppm[0])) for bounds in regions_ppm)

    p0, p1 = settings.p0_deg, settings.p1_deg
    from_pivot = (pivot - ppm) / width_ppm
    peak_index = None
    phase_fit = None
    if settings.autophase_method == "lineshape":
        edges = tuple(edge for bounds in regions for edge in bounds) or ppm_limits
        p0, p1, phase_fit = _lineshape_phase(phase_data, ppm, pivot, width_ppm, edges)
    elif settings.autophase_method == "acme":
        if np.all(phase_data == phase_data[0]):
            raise ValueError("Autophase requires a nonconstant spectrum")
        magnitude = np.abs(phase_data)
        peak_index = np.argmax(magnitude)
        # ACME's score depends on intensity. Normalize only its input copy.
        scaled = phase_data / magnitude[peak_index]
        # nmrglue uses the leftmost sample as the pivot and index / size as the ramp.
        _, (left_phase, p1) = autops(
            scaled, settings.autophase_method, p0=p0 + p1 * from_pivot[0], p1=p1,
            return_phases=True, disp=False,
        )
        p0 = left_phase - p1 * from_pivot[0]
    if phase_data is data or peak_index is not None:
        if p0 or p1:
            phase_data *= np.exp(1j * np.deg2rad(p0 + p1 * from_pivot))
        if peak_index is not None and phase_data[peak_index].real < 0:
            p0 += 180.0
            phase_data *= -1
    if phase_data is not data and (p0 or p1):
        data *= np.exp(1j * np.deg2rad(p0 + p1 * from_pivot))

    windows = tuple(tuple(window) for window in settings.baseline_windows_ppm)
    applied = replace(settings, p0_deg=float(p0), p1_deg=float(p1), pivot_ppm=pivot, baseline_windows_ppm=windows)
    spectrum = Spectrum(
        data=data, hz=hz, ppm=ppm, source=raw.source, intensity_note=raw.intensity_note,
        delay=correction.decision, delay_removal_method=correction.method, tail_trim_points=correction.tail_trim_points,
        acquired_points=raw.complex_points, transformed_points=size, dwell_s=dwell, spectral_width_hz=sw,
        observe_frequency_mhz=observe, reference_ppm=reference, ppm_limits=ppm_limits, settings=applied,
        baseline_coefficients=None, first_point_weight=first_point_weight,
        transform="nmrglue fft_positive with negated frequency bins", lineshape_phase_fit=phase_fit,
        baseline_limit=baseline_limit,
    )
    if settings.smooth_baseline:
        return _apply_automatic_baseline(spectrum, regions, baseline_shape=uses_quantitative_settings(settings))
    return _apply_linear_baseline(spectrum, regions) if windows else spectrum


def baseline_shape_areas(spectrum: Spectrum, bounds_ppm: tuple[float, float]) -> BaselineShapeAreas | None:
    """Look up candidate areas recorded for a declared region under the quantitative profile."""
    if spectrum.smooth_baseline_fit is not None:
        bounds = tuple(sorted(bounds_ppm))
        return next((region.shape_areas for region in spectrum.smooth_baseline_fit.regions
                     if region.bounds_ppm == bounds), None)
    return None


def integrate(spectrum: Spectrum, bounds_ppm: tuple[float, float]) -> Integral:
    """Integrate signed real intensity in increasing Hz within the computed grid."""
    low, high = sorted(bounds_ppm)
    if not low < high:
        raise ValueError("ppm bounds must span a nonzero interval")
    hz = spectrum.hz[::-1]
    # Convert the affine ppm axis to Hz using exact grid endpoints.
    grid_ppm = (spectrum.ppm[-1], spectrum.ppm[0])
    grid_hz = (hz[0], hz[-1])
    low_hz, high_hz = np.interp((low, high), grid_ppm, grid_hz, left=-np.inf, right=np.inf)
    if not hz[0] <= low_hz <= high_hz <= hz[-1]:
        raise ValueError("Integration region reaches the band edge of the computed grid")
    bounds_hz = (low_hz, high_hz)
    intensity = spectrum.data.real[::-1]
    endpoint_values = np.interp(bounds_hz, hz, intensity)
    start = hz.searchsorted(bounds_hz[0], side="right")
    stop = hz.searchsorted(bounds_hz[1])
    x = np.concatenate((bounds_hz[:1], hz[start:stop], bounds_hz[1:]))
    y = np.concatenate((endpoint_values[:1], intensity[start:stop], endpoint_values[1:]))
    area = 0.5 * np.dot(np.diff(x), y[:-1] + y[1:])
    return Integral((low, high), bounds_hz, float(area))
