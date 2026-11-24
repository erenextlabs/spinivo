"""Assess stability of a declared integral on the processed spectrum."""

from __future__ import annotations

from dataclasses import dataclass, replace
from cmath import exp
from math import inf, isclose, isfinite, radians
from numbers import Real
from typing import TYPE_CHECKING

from .acquisition import CheckResult, SignalEvidence, ThresholdProfile, check_snr

if TYPE_CHECKING:
    from .processing import Spectrum
    from .readers.model import RawFID


_REFIT_THRESHOLD_VERSION = "pilot-2"


@dataclass(frozen=True)
class IntegralThresholds:
    version: str = _REFIT_THRESHOLD_VERSION
    max_relative_change: float = 0.01
    min_snr: float = 150.0
    min_noise_points: int = 32
    unexpected_peak_sd: float = 5.0
    max_pair_capture_mismatch: float = 0.005

    def __post_init__(self):
        for field in ("max_relative_change", "min_snr", "unexpected_peak_sd", "max_pair_capture_mismatch"):
            value = getattr(self, field)
            if not isinstance(value, Real) or isinstance(value, bool) or not 0 < value < inf:
                raise ValueError(f"{field} must be finite and positive")
        if type(self.min_noise_points) is not int or self.min_noise_points <= 0:
            raise ValueError("min_noise_points must be a positive integer")
        if not isinstance(self.version, str) or not self.version:
            raise ValueError("version must be a nonempty string")


PILOT = IntegralThresholds()
# Rounded upper Tukey fence of the July FID02 regions, pending independent validation.
_BASELINE_SENSITIVITY_LIMIT = 0.017
_BASELINE_SENSITIVITY_FLAG = "Baseline position or placement sensitivity exceeds the limit"


def _spike_anchor_lock_applied(spectrum: Spectrum) -> bool:
    return spectrum.smooth_baseline_fit is not None and spectrum.settings.lock_spike_anchors


def _region_thresholds(thresholds, lock_spike_anchors):
    if lock_spike_anchors and thresholds.version == _REFIT_THRESHOLD_VERSION:
        return replace(thresholds, version="pilot-3")
    return thresholds


@dataclass(frozen=True)
class BaselineRecovery:
    bounds_ppm: tuple[float, float]
    value: float | None  # Midpoint recovered area divided by added area.
    reason: str | None = None
    limit: float = _BASELINE_SENSITIVITY_LIMIT
    flag: str | None = None
    recovery_25: float | None = None
    recovery_75: float | None = None
    position_sensitivity: float | None = None
    placement_sensitivity: float | None = None


@dataclass(frozen=True)
class SpikeAnchorFlips:
    """Quiet blocks an unlocked spike trial would add or remove, by center ppm."""

    position_percent: int
    added_ppm: tuple[float, ...]
    removed_ppm: tuple[float, ...]


@dataclass(frozen=True)
class LockedBaselineRecovery(BaselineRecovery):
    anchor_flips: tuple[SpikeAnchorFlips, ...] = ()


def check_baseline_recovery(spectrum: Spectrum, bounds_ppm: tuple[float, float]) -> BaselineRecovery:
    """Return the recovery check recorded while processing a declared region."""
    bounds = tuple(sorted(bounds_ppm))
    if spectrum.smooth_baseline_fit is not None:
        for region in spectrum.smooth_baseline_fit.regions:
            if region.bounds_ppm == bounds and region.recovery is not None:
                return region.recovery
    for recovery in spectrum.linear_baseline_recoveries:
        if recovery.bounds_ppm == bounds:
            return recovery
    if spectrum.smooth_baseline_fit is None and spectrum.baseline_coefficients is None:
        return BaselineRecovery(bounds, None, "Spectrum has no baseline correction", limit=spectrum.baseline_limit)
    return BaselineRecovery(bounds, None, "No recovery recorded for these bounds during processing",
                            limit=spectrum.baseline_limit)


def _recovery_trial(spectrum, bounds, center_ppm, nominal_area, refit, position_percent):
    from .processing import integrate

    center_hz = (center_ppm - spectrum.reference_ppm) * spectrum.observe_frequency_mhz
    # A 2 Hz FWHM absorptive line with 1 Hz half-width, sampled on the phased grid.
    line = 1 / (1 + (spectrum.hz - center_hz) ** 2)
    sampled_area = integrate(replace(spectrum, data=line), bounds).area
    added_area = 0.10 * nominal_area
    if not isfinite(sampled_area) or sampled_area <= 0 or added_area <= 0:
        raise ValueError("Cannot scale the addition to a positive in-window area")
    trial = replace(spectrum, data=spectrum.data + line * (added_area / sampled_area))
    area = integrate(refit(trial, position_percent), bounds).area
    value = (area - nominal_area) / added_area
    if not isfinite(value):
        raise ValueError("Recovered area ratio is nonfinite")
    return value


def _measure_baseline_recovery(spectrum, bounds, nominal_area, refit, straight_area, *,
                               limit=_BASELINE_SENSITIVITY_LIMIT):
    if not isfinite(nominal_area) or nominal_area <= 0:
        return BaselineRecovery(bounds, None, "Nominal area must be finite and positive", limit=limit)
    low, high = bounds
    midpoint = sum(bounds) / 2
    width = high - low
    reasons = []
    values = []
    for position, label, center in ((50, "Midpoint", midpoint), (25, "25%", low + 0.25 * width),
                                    (75, "75%", low + 0.75 * width)):
        try:
            values.append(_recovery_trial(spectrum, bounds, center, nominal_area, refit, position))
        except ValueError as error:
            reasons.append(f"{label} recovery: {error}")
            values.append(None)
    midpoint_recovery, recovery_25, recovery_75 = values
    position = None
    if recovery_25 is not None and recovery_75 is not None:
        position = abs(recovery_25 - recovery_75)
        if not isfinite(position):
            position = None
            reasons.append("Position sensitivity is nonfinite")
    placement = None
    if straight_area is None:
        reasons.append("Missing baseline anchor bracket for straight-line placement")
    else:
        placement = abs(straight_area - nominal_area) / nominal_area
        if not isfinite(placement):
            placement = None
            reasons.append("Placement sensitivity is nonfinite")
    exceeds_limit = any(score is not None and score > limit for score in (position, placement))
    return BaselineRecovery(bounds, midpoint_recovery, ". ".join(reasons) or None, limit=limit,
                            flag="baseline_sensitive" if exceeds_limit else None,
                            recovery_25=recovery_25, recovery_75=recovery_75,
                            position_sensitivity=position, placement_sensitivity=placement)


@dataclass(frozen=True)
class BaselineDetails:
    anchor_distances_ppm: tuple[float | None, float | None]
    recovery: BaselineRecovery


@dataclass(frozen=True)
class BoundsTrial:
    move: str
    bounds_ppm: tuple[float, float]
    relative_change: float | None
    reason: str | None


@dataclass(frozen=True)
class PhaseTrial:
    peak_ppm: float
    residual_deg: float
    relative_change: float | None
    reason: str | None


@dataclass(frozen=True)
class IntegralObservation:
    check: str
    status: str
    value: float | None
    unit: str
    limit: float | None
    reason: str
    threshold_version: str
    details: BaselineDetails | tuple[BoundsTrial, ...] | tuple[PhaseTrial, ...] | CheckResult | None = None


@dataclass(frozen=True)
class RegionStability:
    bounds_ppm: tuple[float, float]
    nominal_area: float | None
    area: IntegralObservation
    baseline: IntegralObservation
    bounds: IntegralObservation
    phase: IntegralObservation
    snr: IntegralObservation
    thresholds: IntegralThresholds


def _relative_flags(value, limit, flag=None):
    flags = ("Relative area change exceeds the limit",) if value is not None and value > limit else ()
    return (*flags, flag) if flag else flags


def _relative_observation(check, value, thresholds, details=None, *, reason=None, flag=None):
    limit = thresholds.max_relative_change
    if flags := _relative_flags(value, limit, flag):
        status = "flag"
        reason = ". ".join((*flags, reason) if reason else flags)
    elif reason:
        status = "cannot be assessed"
    else:
        status, reason = "pass", "Relative area change meets the limit"
    return IntegralObservation(check, status, value, "fraction", limit, reason, thresholds.version, details)


def _baseline_observation(spectrum, bounds, area, area_reason, thresholds):
    fit = spectrum.smooth_baseline_fit
    region = next((r for r in fit.regions if r.bounds_ppm == bounds), None) if fit is not None else None
    recovery = check_baseline_recovery(spectrum, bounds)
    distances = region.anchor_distances_ppm if region is not None else (None, None)
    gaps = [f"Missing {side} baseline anchor" for side, d in zip(("left", "right"), distances)
            if d is None or not isfinite(d)]
    alternate = region.half_block_area if region is not None else None
    if alternate is None or not isfinite(alternate):
        gaps.append("No finite half-block area for these bounds")
    if area_reason:
        gaps.append(area_reason)
    value = None
    if not gaps:
        value = abs(alternate - area) / abs(area)
        if not isfinite(value):
            value = None
            gaps.append("Half-block relative area change is nonfinite")
    if recovery.reason:
        gaps.append(recovery.reason)
    flag = _BASELINE_SENSITIVITY_FLAG if recovery.flag else None
    return _relative_observation("baseline", value, thresholds, BaselineDetails(distances, recovery),
                                 reason=". ".join(gaps) or None, flag=flag)


def _trial_change(spectrum, bounds, area):
    from .processing import integrate

    try:
        change = (integrate(spectrum, bounds).area - area) / abs(area)
    except ValueError as error:
        return None, str(error)
    return (change, None) if isfinite(change) else (None, "Relative area change is nonfinite")


def _trial_observation(check, trials, thresholds):
    value = max((abs(t.relative_change) for t in trials if t.relative_change is not None), default=None)
    gaps = dict.fromkeys(t.reason for t in trials if t.reason)
    reason = f"Incomplete {check} trials. " + ". ".join(gaps) if gaps else None
    return _relative_observation(check, value, thresholds, tuple(trials), reason=reason)


def _bounds_observation(spectrum, bounds, area, area_reason, thresholds):
    low, high = bounds
    step = max(spectrum.spectral_width_hz / spectrum.acquired_points / spectrum.observe_frequency_mhz,
               0.05 * (high - low))
    trials = []
    for move, trial in (("low_out", (low - step, high)), ("low_in", (low + step, high)),
                        ("high_in", (low, high - step)), ("high_out", (low, high + step))):
        start, end = trial
        change, reason = None, area_reason
        if start >= end:
            reason = "Trial bounds collapse"
        elif not spectrum.ppm_limits[0] <= start < end <= spectrum.ppm_limits[1]:
            reason = "Trial bounds reach outside the spectral width"
        elif not area_reason:
            change, reason = _trial_change(spectrum, trial, area)
        trials.append(BoundsTrial(move, trial, change, reason))
    if area_reason:
        return _relative_observation("bounds", None, thresholds, tuple(trials), reason=area_reason)
    return _trial_observation("bounds", trials, thresholds)


def _phase_observation(spectrum, bounds, area, area_reason, thresholds):
    fit = spectrum.lineshape_phase_fit
    reason = area_reason
    if not reason:
        if fit is None:
            reason = "No lineshape phase fit"
        elif fit.method == "single_peak":
            reason = "The single_peak fit sets the residual to zero by construction"
        else:
            peaks = [peak for peak in fit.peaks if bounds[0] <= peak.ppm <= bounds[1]]
            if not peaks:
                reason = "No accepted phase peak inside the bounds"
    if reason:
        return _relative_observation("phase", None, thresholds, reason=reason)

    baseline = spectrum.smooth_baseline_curve
    if baseline is None:
        if spectrum.smooth_baseline_fit is not None or spectrum.settings.smooth_baseline:
            return _relative_observation("phase", None, thresholds, reason="Subtracted smooth baseline is missing")
        baseline = 0.0
        if spectrum.baseline_coefficients is not None:
            intercept, slope = spectrum.baseline_coefficients
            baseline = intercept + slope * spectrum.ppm
    uncorrected = spectrum.data + baseline
    settings = spectrum.settings
    pivot = settings.pivot_ppm if settings.pivot_ppm is not None else spectrum.reference_ppm
    width = spectrum.spectral_width_hz / spectrum.observe_frequency_mhz
    trials = []
    for peak in peaks:
        # A transferred fit supplies residuals from the phase-estimation trace.
        # The area trial uses the quantitative trace and its own baseline.
        applied = settings.p0_deg + settings.p1_deg * (pivot - peak.ppm) / width
        residual = (peak.phase_deg - applied + 180) % 360 - 180
        if isfinite(residual):
            rotated = replace(spectrum, data=uncorrected * exp(1j * radians(residual)) - baseline)
            change, reason = _trial_change(rotated, bounds, area)
        else:
            change, reason = None, "Phase residual is nonfinite"
        trials.append(PhaseTrial(peak.ppm, residual, change, reason))
    return _trial_observation("phase", trials, thresholds)


def _snr_input_reason(spectrum, raw, noise_bounds):
    if raw is None or noise_bounds is None:
        return "Missing raw acquisition" if raw is None else "Missing noise window"
    if raw.complex_points != spectrum.acquired_points:
        return "Raw complex_points does not match the processed spectrum"
    for field in ("spectral_width_hz", "observe_frequency_mhz"):
        value = getattr(raw, field)
        if value is None or not isclose(value, getattr(spectrum, field), rel_tol=1e-12):
            return f"Raw {field} does not match the processed spectrum"
    if raw.source and spectrum.source and raw.source != spectrum.source:
        return "Raw source does not match the processed spectrum"
    return None


def _snr_observation(spectrum, bounds, raw, noise_bounds, thresholds):
    from .measurement import measure_region

    limit, version = thresholds.min_snr, thresholds.version
    if reason := _snr_input_reason(spectrum, raw, noise_bounds):
        return IntegralObservation("snr", "cannot be assessed", None, "ratio", limit, reason, version)
    measured = measure_region(spectrum, "integral region", bounds, noise_bounds)
    signal = SignalEvidence(name=measured.name, measurement=measured.to_signal_measurement())
    result = check_snr(raw, signal, ThresholdProfile(min_snr=limit, min_noise_points=thresholds.min_noise_points))
    return IntegralObservation("snr", result.status, result.value, result.unit, result.limit,
                               result.reason, version, result)


def check_region_stability(spectrum: Spectrum, bounds: tuple[float, float], *, raw: RawFID | None = None,
                           noise_bounds_ppm: tuple[float, float] | None = None,
                           thresholds: IntegralThresholds = PILOT) -> RegionStability:
    """Return independent observations using the nominal processed trace and signed area."""
    from .processing import integrate

    thresholds = _region_thresholds(thresholds, _spike_anchor_lock_applied(spectrum))
    bounds = tuple(sorted(bounds))
    area, reason = None, None
    try:
        area = integrate(spectrum, bounds).area
        if not isfinite(area):
            reason = "Nominal area is nonfinite"
        elif area <= 0:
            reason = "Nominal area is nonpositive"
    except ValueError as error:
        reason = str(error)
    nominal = IntegralObservation("area", "cannot be assessed" if reason else "pass", area, "intensity*Hz",
                                  None, reason or "Nominal area is finite and positive", thresholds.version)
    return RegionStability(
        bounds, area, nominal, _baseline_observation(spectrum, bounds, area, reason, thresholds),
        _bounds_observation(spectrum, bounds, area, reason, thresholds),
        _phase_observation(spectrum, bounds, area, reason, thresholds),
        _snr_observation(spectrum, bounds, raw, noise_bounds_ppm, thresholds),
        thresholds,
    )
