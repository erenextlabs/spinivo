"""Check assigned Lorentzian components and combine integral observations."""

from __future__ import annotations

from dataclasses import dataclass
from math import atan, inf, isfinite, nan, pi
from typing import Literal

import numpy as np

from .integral_checks import (
    PILOT, BaselineDetails, IntegralObservation, IntegralThresholds, RegionStability,
    _BASELINE_SENSITIVITY_FLAG, _region_thresholds, _relative_flags, _spike_anchor_lock_applied, check_region_stability,
)
from .measurement import measure_region
from .processing import Spectrum
from .readers.model import RawFID


@dataclass(frozen=True)
class Line:
    """Whole-line area in intensity*Hz and effective FWHM on the processed trace."""

    ppm: float
    area: float
    fwhm_hz: float

    def __post_init__(self):
        if not all(map(isfinite, (self.ppm, self.area, self.fwhm_hz))) or min(self.area, self.fwhm_hz) <= 0:
            raise ValueError("Line requires finite ppm, positive area and positive FWHM")


@dataclass(frozen=True)
class ComponentEvidence:
    target_lines: tuple[Line, ...] = ()
    target_source: str | None = None
    foreign_lines: tuple[Line, ...] = ()
    foreign_source: str | None = None
    satellite_convention: Literal["include", "exclude", "absent"] | None = None
    coupling_hz: float | None = None
    satellite_fraction: float | None = 0.0055
    satellite_source: str | None = None

    def __post_init__(self):
        if self.satellite_convention not in (None, "include", "exclude", "absent"):
            raise ValueError("Unknown satellite convention")
        if self.coupling_hz is not None and not 0 < self.coupling_hz < inf:
            raise ValueError("coupling_hz must be finite and positive")
        if self.satellite_fraction is not None and not 0 <= self.satellite_fraction <= 0.5:
            raise ValueError("satellite_fraction must be between zero and one half")


@dataclass(frozen=True)
class SatelliteBand:
    ppm: float
    fwhm_hz: float
    captured: float
    status: str


@dataclass(frozen=True)
class ComponentObservation:
    check: str
    value: float | None
    limit: float
    threshold_version: str
    flags: tuple[str, ...] = ()
    gaps: tuple[str, ...] = ()
    bands: tuple[SatelliteBand, ...] = ()
    unit: str = "fraction"

    @property
    def status(self) -> str:
        return _status(self.flags, self.gaps)

    @property
    def reason(self) -> str:
        return ". ".join((*self.flags, *self.gaps)) or "Meets the pilot check"


@dataclass(frozen=True)
class RegionResult(RegionStability):
    """Region observations with the source of the processed spectrum."""

    evidence: ComponentEvidence
    neighbour: ComponentObservation
    capture: ComponentObservation
    satellites: ComponentObservation
    status: str
    flags: tuple[tuple[str, str], ...]
    gaps: tuple[tuple[str, str], ...]
    requested_thresholds: IntegralThresholds
    lock_spike_anchors: bool = False
    run: str | None = None

    @property
    def observations(self) -> tuple[IntegralObservation | ComponentObservation, ...]:
        return (self.area, self.baseline, self.bounds, self.phase, self.snr,
                self.neighbour, self.capture, self.satellites)


def _status(flags, gaps):
    if flags:
        return "flag"
    return "cannot be assessed" if gaps else "pass"


def _captured(center, width, bounds, scale):
    if not isfinite(center) or not 0 < width < inf:
        return nan
    low, high = bounds
    return (atan(2 * scale * (high - center) / width) - atan(2 * scale * (low - center) / width)) / pi


def _target_weights(lines):
    largest = max(line.area for line in lines)
    return [line.area / largest for line in lines]


def _target_gaps(evidence):
    values = (evidence.target_lines, (evidence.target_source or "").strip())
    return tuple(f"Missing target {name}" for name, value in zip(("lines", "source"), values) if not value)


def _capture(evidence, scale, bounds, bounds_gap, thresholds):
    gaps = _target_gaps(evidence) + bounds_gap
    limit = thresholds.max_relative_change
    if gaps:
        return ComponentObservation("capture", None, limit, thresholds.version, gaps=gaps)
    lines = evidence.target_lines
    fractions = [_captured(line.ppm, line.fwhm_hz, bounds, scale) for line in lines]
    loss = 1 - float(np.average(fractions, weights=_target_weights(lines)))
    if not isfinite(loss):
        return ComponentObservation("capture", None, limit, thresholds.version, gaps=("Capture loss is nonfinite",))
    flags = ("Target capture loss exceeds the limit",) if loss > limit else ()
    return ComponentObservation("capture", loss, limit, thresholds.version, flags)


def _satellites(evidence, scale, bounds, bounds_gap, thresholds):
    limit = thresholds.max_relative_change
    convention = evidence.satellite_convention
    fraction = evidence.satellite_fraction
    gaps = tuple(reason for missing, reason in (
        (convention is None, "Missing satellite convention"),
        (not (evidence.satellite_source or "").strip(), "Missing satellite source"),
        (convention in ("include", "exclude") and evidence.coupling_hz is None, "Missing satellite coupling"),
    ) if missing)
    if convention in ("include", "exclude"):
        gaps += _target_gaps(evidence) + bounds_gap
    fraction_gap = ("Missing satellite fraction",) if convention in ("include", "exclude") and fraction is None else ()
    if gaps or convention == "absent":
        return ComponentObservation("satellites", None if gaps else 0.0, limit, thresholds.version,
                                    gaps=gaps + fraction_gap)

    lines = evidence.target_lines
    offset = evidence.coupling_hz / (2 * scale)
    bands, errors = [], []
    asymmetric = False
    for line in lines:
        pair = []
        for ppm in (line.ppm - offset, line.ppm + offset):
            captured = _captured(ppm, line.fwhm_hz, bounds, scale)
            status = "inside" if captured >= 0.95 else "outside" if captured <= 0.05 else "clipped"
            pair.append(SatelliteBand(ppm, line.fwhm_hz, captured, status))
        bands.extend(pair)
        errors.append(sum(1 - band.captured if convention == "include" else band.captured for band in pair))
        if pair[0].status != pair[1].status:
            asymmetric = True
    if fraction_gap:
        return ComponentObservation("satellites", None, limit, thresholds.version, gaps=fraction_gap, bands=tuple(bands))
    error = fraction * float(np.average(errors, weights=_target_weights(lines)))
    if not isfinite(error):
        return ComponentObservation("satellites", None, limit, thresholds.version,
                                    gaps=("Satellite area error is nonfinite",))
    flags = ("Satellite area error exceeds the limit",) if error > limit else ()
    if asymmetric:
        flags += ("satellite_asymmetric",)
    return ComponentObservation("satellites", error, limit, thresholds.version, flags, bands=tuple(bands))


def _resolved_peak(real, index, cutoff):
    # The rubric screens resolved maxima, so noise ripples on assigned wings must not count.
    # Require prominence on both flanks.
    height = real[index]
    for flank in (real[index - 1::-1], real[index + 1:]):
        for point in flank:
            if point < height - cutoff:
                break
            if point > height:
                return False
        else:
            return False
    return True


def _unexpected_peaks(spectrum, bounds, noise_bounds, assigned, thresholds):
    if noise_bounds is None:
        return ("Missing noise window",)
    measured = measure_region(spectrum, "integral region", bounds, noise_bounds)
    noise = measured.noise_sd
    if noise is None or not isfinite(noise):
        return (measured.noise_sd_reason or "Noise SD is nonfinite",)
    ppm, real = spectrum.ppm[1:-1], spectrum.data.real
    cutoff = thresholds.unexpected_peak_sd * noise
    middle = real[1:-1]
    peaks = np.flatnonzero((bounds[0] <= ppm) & (ppm <= bounds[1]) & (middle > cutoff)
                           & (middle > real[:-2]) & (middle >= real[2:]))
    scale = spectrum.observe_frequency_mhz
    return tuple(f"Unexpected peak at {ppm[i]:.5f} ppm" for i in peaks
                 if all(abs(ppm[i] - line.ppm) * scale > line.fwhm_hz for line in assigned)
                 and _resolved_peak(real, i + 1, cutoff))


def _neighbour(spectrum, stability, evidence, bands, noise_bounds, thresholds):
    bounds, area = stability.bounds_ppm, stability.nominal_area
    limit = thresholds.max_relative_change
    gaps, flags, q = [], (), None
    if evidence.foreign_lines and not (evidence.foreign_source or "").strip():
        gaps.append("Missing foreign source")
    if stability.area.status != "pass":
        gaps.append(stability.area.reason)
    if not gaps:
        scale = spectrum.observe_frequency_mhz
        q = sum(line.area * _captured(line.ppm, line.fwhm_hz, bounds, scale)
                for line in evidence.foreign_lines) / abs(area)
        if q > limit:
            flags = ("Foreign area fraction exceeds the limit",)
    gaps.extend(_unexpected_peaks(spectrum, bounds, noise_bounds,
                                  (*evidence.target_lines, *evidence.foreign_lines, *bands), thresholds))
    return ComponentObservation("neighbour", q, limit, thresholds.version, flags, tuple(gaps))


def _stability_causes(observation):
    if observation.status != "flag":
        return (), (observation.reason,) if observation.status == "cannot be assessed" else ()
    details = observation.details
    sensitivity = _BASELINE_SENSITIVITY_FLAG if isinstance(details, BaselineDetails) and details.recovery.flag else None
    flags = _relative_flags(observation.value, observation.limit, sensitivity)
    if not flags:
        return (observation.reason,), ()  # S/N uses its own check.
    reason = observation.reason.removeprefix(". ".join(flags)).removeprefix(". ")
    return flags, (reason,) if reason else ()


def check_region(spectrum: Spectrum, bounds: tuple[float, float], evidence: ComponentEvidence, *,
                 raw: RawFID | None = None, noise_bounds_ppm: tuple[float, float] | None = None,
                 thresholds: IntegralThresholds = PILOT) -> RegionResult:
    """Check all eight observations and retain supported flags and evidence gaps."""
    stability = check_region_stability(spectrum, bounds, raw=raw, noise_bounds_ppm=noise_bounds_ppm,
                                       thresholds=thresholds)
    lock_spike_anchors = _spike_anchor_lock_applied(spectrum)
    bounds, scale = stability.bounds_ppm, spectrum.observe_frequency_mhz
    bounds_gap = (stability.area.reason,) if stability.nominal_area is None else ()
    capture = _capture(evidence, scale, bounds, bounds_gap, stability.thresholds)
    satellites = _satellites(evidence, scale, bounds, bounds_gap, stability.thresholds)
    neighbour = _neighbour(spectrum, stability, evidence, satellites.bands, noise_bounds_ppm, stability.thresholds)
    stability_observations = (stability.area, stability.baseline, stability.bounds, stability.phase, stability.snr)
    causes = [(obs.check, *_stability_causes(obs)) for obs in stability_observations]
    causes.extend((obs.check, obs.flags, obs.gaps) for obs in (neighbour, capture, satellites))
    flags = tuple((check, reason) for check, reasons, _ in causes for reason in reasons)
    gaps = tuple((check, reason) for check, _, reasons in causes for reason in reasons)
    return RegionResult(**vars(stability), evidence=evidence, neighbour=neighbour, capture=capture,
                        satellites=satellites, status=_status(flags, gaps), flags=flags, gaps=gaps,
                        requested_thresholds=thresholds, lock_spike_anchors=lock_spike_anchors,
                        run=None if spectrum.source is None else str(spectrum.source))


def check_capture_pair(analyte: RegionResult, standard: RegionResult,
                       thresholds: IntegralThresholds = PILOT) -> ComponentObservation:
    """Check the relative analyte-to-standard response from the recorded target losses."""
    limit = thresholds.max_pair_capture_mismatch
    if analyte.lock_spike_anchors != standard.lock_spike_anchors:
        return ComponentObservation("capture_pair", None, limit, thresholds.version,
                                    gaps=("Analyte and standard use different spike anchor lock modes",))
    thresholds = _region_thresholds(thresholds, analyte.lock_spike_anchors)
    analyte_loss, standard_loss = analyte.capture.value, standard.capture.value
    gaps = tuple(f"Missing or unusable {side} capture"
                 for side, loss in (("analyte", analyte_loss), ("standard", standard_loss))
                 if loss is None or not 0 <= loss <= 1 or (side == "standard" and loss == 1))
    if gaps:
        return ComponentObservation("capture_pair", None, limit, thresholds.version, gaps=gaps)
    ratio = (1 - analyte_loss) / (1 - standard_loss)
    flags = () if 1 - limit <= ratio <= 1 + limit else ("Pair capture mismatch exceeds the limit",)
    return ComponentObservation("capture_pair", abs(ratio - 1), limit, thresholds.version, flags)
