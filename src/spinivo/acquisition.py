"""Assess qNMR acquisition settings and measured signal quality."""

from collections.abc import Iterable
from dataclasses import dataclass
from math import exp, expm1, isclose, isfinite, radians, sin
from numbers import Real

from .readers.model import RawFID


@dataclass(frozen=True)
class ThresholdProfile:
    """Pilot proposals for the acquisition checks."""

    max_recovery_loss: float = 0.01
    min_snr: float = 150.0
    min_points_per_fwhm: float = 5.0
    min_noise_points: int = 32

    def __post_init__(self):
        if type(self.min_noise_points) is not int or self.min_noise_points <= 0:
            raise ValueError("min_noise_points must be a positive integer")


PILOT_THRESHOLDS = ThresholdProfile()


@dataclass(frozen=True)
class RecoveryEvidence:
    """Applicable pulse calibration with its traceable source and optional extra recovery for one run."""

    calibrated_flip_angle_deg: float | None = None
    calibrated_90_us: float | None = None
    calibration_source: str | None = None
    extra_recovery_s: float | None = 0.0
    extra_recovery_source: str | None = None


@dataclass(frozen=True)
class SignalMeasurement:
    """Qualified height, sample noise SD, noise window and linewidth."""

    height: float | None = None
    noise_sd: float | None = None
    fwhm_hz: float | None = None
    source: str = "signal measurement"
    noise_bounds_ppm: tuple[float, float] | None = None
    transformed_points: int | None = None
    fwhm_reason: str | None = None


@dataclass(frozen=True)
class SignalEvidence:
    name: str | None
    t1_s: float | None = None
    t1_source: str | None = None
    measurement: SignalMeasurement = SignalMeasurement()


@dataclass(frozen=True)
class CheckInput:
    name: str
    value: float | int | str | tuple[float, float] | None
    source: str


@dataclass(frozen=True)
class CheckResult:
    run: str
    signal: str | None
    check: str
    status: str
    value: float | None
    limit: float
    unit: str
    reason: str
    inputs: tuple[CheckInput, ...]


@dataclass(frozen=True)
class PairAcquisitionResults:
    sample: tuple[CheckResult, ...]
    calibrant: tuple[CheckResult, ...]


_PULSE_PROGRAMS = frozenset(("zg", "zg30", "zgig", "zgig30", "single_pulse", "single_pulse_dec"))


def _has_text(value):
    return isinstance(value, str) and bool(value.strip())


def _raw_input(raw, name):
    return CheckInput(name, getattr(raw, name), raw.source_fields.get(name, f"RawFID.{name}"))


def _positive_inputs(*inputs):
    for item in inputs:
        if item.value is None:
            return f"{item.name} is missing"
        if not isinstance(item.value, Real) or not isfinite(item.value) or item.value <= 0:
            return f"{item.name} must be finite and positive"
    return None


def _result(raw, signal, check, limit, unit, inputs, *, value=None, reason=None):
    if value is None:
        status = "cannot be assessed"
    else:
        passed = value <= limit if check == "recovery" else value >= limit
        status = "pass" if passed else "flag"
        reason = "Meets the pilot limit" if passed else "Outside the pilot limit"
        reason += " for the steady-state bound" if check == "recovery" else ""
    return CheckResult(str(raw.source), signal.name, check, status, value, limit, unit, reason, inputs)


def check_recovery(raw: RawFID, signal: SignalEvidence, evidence: RecoveryEvidence = RecoveryEvidence(),
                   thresholds: ThresholdProfile = PILOT_THRESHOLDS) -> CheckResult:
    """Bound recovery loss at steady state for a pulse train starting at equilibrium."""
    program = _raw_input(raw, "pulse_program")
    d1 = _raw_input(raw, "relaxation_delay_s")
    aq = _raw_input(raw, "acquisition_time_s")
    width = _raw_input(raw, "pulse_width_us")
    t1 = CheckInput("t1_s", signal.t1_s, signal.t1_source)
    extra = evidence.extra_recovery_s
    extra_source = evidence.extra_recovery_source or ("default" if extra == 0 else "assay evidence")
    inputs = (program, d1, aq, width, t1, CheckInput("extra_recovery_s", extra, extra_source),
              _raw_input(raw, "scans"), _raw_input(raw, "dummy_scans"))
    angle = evidence.calibrated_flip_angle_deg
    cal90 = CheckInput("calibrated_90_us", evidence.calibrated_90_us, evidence.calibration_source)
    angle_input = CheckInput("calibrated_flip_angle_deg", angle, evidence.calibration_source)
    inputs += (angle_input,) if angle is not None else (cal90,)
    reason = _positive_inputs(d1, aq, width, t1)
    if reason is None and not _has_text(signal.t1_source):
        reason = "t1_source is missing"
    if program.value is None or not str(program.value).strip():
        reason = "pulse_program is missing, so the excitation scheme is unknown"
    elif program.value not in _PULSE_PROGRAMS:
        reason = f"pulse_program {program.value!r} is not a supported 1D pulse-acquire program"
    elif not reason:
        if extra is None or not isfinite(extra) or extra < 0:
            reason = "extra_recovery_s must be finite and nonnegative"
        elif angle is None:
            reason = _positive_inputs(cal90)
            if reason is None:
                if raw.vendor == "bruker" and program.value in ("zg30", "zgig30"):
                    width = CheckInput("effective_pulse_width_us", width.value * 0.33,
                                       "Bruker zg30/zgig30 pulse program: p1*0.33")
                    inputs += (width,)
                angle = 90 * (width.value / cal90.value)
                inputs += (CheckInput("calibrated_flip_angle_deg", angle,
                                      f"90 * {width.name} / calibrated_90_us"),)
            elif cal90.value is None:
                reason = "calibrated_flip_angle_deg or calibrated_90_us is required"
        if reason is None and not _has_text(evidence.calibration_source):
            reason = "calibration_source is missing"
        if reason is None and not 0 < angle <= 180:
            reason = "calibrated_flip_angle_deg must be in (0, 180]"

    loss = None
    if reason is None:
        tr = d1.value + aq.value + extra
        inputs += (CheckInput("tr_s", tr, "relaxation_delay_s + acquisition_time_s + extra_recovery_s"),)
        try:
            scaled_tr = tr / t1.value
            finite_ratio = isfinite(scaled_tr)
        except (ArithmeticError, RuntimeWarning):
            finite_ratio = False
        if not finite_ratio:
            reason = "tr_s / t1_s must be finite to evaluate recovery"
        else:
            # Evaluate 1 - R without subtracting two nearly equal responses.
            saturation = exp(-scaled_tr) * (2 * sin(radians(angle) / 2) ** 2)
            denominator = saturation - expm1(-scaled_tr)
            if denominator > 0:
                loss = saturation / denominator
            else:
                reason = "calibrated_flip_angle_deg and tr_s / t1_s are too small to evaluate"
    return _result(raw, signal, "recovery", thresholds.max_recovery_loss, "fraction", inputs,
                   value=loss, reason=reason)


def check_snr(raw: RawFID, signal: SignalEvidence, thresholds: ThresholdProfile = PILOT_THRESHOLDS) -> CheckResult:
    """Assess peak height relative to noise measured on the acquired grid."""
    m = signal.measurement
    height = CheckInput("height", m.height, m.source)
    noise = CheckInput("noise_sd", m.noise_sd, m.source)
    bounds = m.noise_bounds_ppm
    observe = _raw_input(raw, "observe_frequency_mhz")
    sw, points = _raw_input(raw, "spectral_width_hz"), _raw_input(raw, "complex_points")
    inputs = (height, noise, CheckInput("noise_bounds_ppm", bounds, m.source), observe, sw, points)
    reason = _positive_inputs(height, noise, observe, sw, points)
    if reason is None and (not isinstance(bounds, tuple) or len(bounds) != 2 or bounds[0] == bounds[1]
                           or not all(isinstance(bound, Real) and isfinite(bound) for bound in bounds)):
        reason = "noise_bounds_ppm requires two distinct finite values"
    if reason is None:
        window = abs(bounds[1] - bounds[0]) * observe.value
        count = window / (sw.value / points.value)
        inputs += (CheckInput("window_hz", window, "noise_bounds_ppm and observe_frequency_mhz"),
                   CheckInput("noise_points", count, "window_hz / (spectral_width_hz / complex_points)"))
        minimum = thresholds.min_noise_points
        if window > sw.value:
            reason = "noise_bounds_ppm window_hz exceeds spectral_width_hz"
        elif count < minimum and not isclose(count, minimum, rel_tol=1e-9):  # ppm endpoint roundoff
            reason = f"noise_points must be at least {minimum} acquired points"
    value = height.value / noise.value if reason is None else None
    return _result(raw, signal, "snr", thresholds.min_snr, "ratio", inputs, value=value, reason=reason)


def check_resolution(raw: RawFID, signal: SignalEvidence,
                     thresholds: ThresholdProfile = PILOT_THRESHOLDS) -> CheckResult:
    """Assess points across an isolated line's zero-broadening FWHM."""
    sw = _raw_input(raw, "spectral_width_hz")
    points = _raw_input(raw, "complex_points")
    m = signal.measurement
    fwhm = CheckInput("fwhm_hz", m.fwhm_hz, m.source)
    inputs = (sw, points, fwhm)
    reason = _positive_inputs(*inputs)
    value = None
    if reason is None:
        spacing = sw.value / points.value
        inputs += (CheckInput("spacing_hz", spacing, "spectral_width_hz / complex_points"),)
        if m.transformed_points is not None and m.transformed_points > points.value:
            inputs += (CheckInput("transformed_points", m.transformed_points, m.source),
                       CheckInput("transformed_spacing_hz", sw.value / m.transformed_points,
                                  "spectral_width_hz / transformed_points"))
        value = fwhm.value / spacing
    elif reason == "fwhm_hz is missing" and m.fwhm_reason:
        reason += f": {m.fwhm_reason}"
    return _result(raw, signal, "resolution", thresholds.min_points_per_fwhm, "points / FWHM", inputs,
                   value=value, reason=reason)


def check_acquisition(raw: RawFID, signals: Iterable[SignalEvidence],
                      evidence: RecoveryEvidence = RecoveryEvidence(),
                      thresholds: ThresholdProfile = PILOT_THRESHOLDS) -> tuple[CheckResult, ...]:
    """Return all three independent checks for each supplied signal."""
    return tuple(result for signal in signals for result in (
        check_recovery(raw, signal, evidence, thresholds),
        check_snr(raw, signal, thresholds),
        check_resolution(raw, signal, thresholds),
    ))


def check_acquisition_pair(sample: RawFID, sample_signals: Iterable[SignalEvidence],
                           calibrant: RawFID, calibrant_signals: Iterable[SignalEvidence], *,
                           sample_evidence: RecoveryEvidence = RecoveryEvidence(),
                           calibrant_evidence: RecoveryEvidence = RecoveryEvidence(),
                           thresholds: ThresholdProfile = PILOT_THRESHOLDS) -> PairAcquisitionResults:
    """Assess sample and calibrant using their own measurements and evidence."""
    return PairAcquisitionResults(
        check_acquisition(sample, sample_signals, sample_evidence, thresholds),
        check_acquisition(calibrant, calibrant_signals, calibrant_evidence, thresholds),
    )
