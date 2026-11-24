"""Assess same-tube or external-calibration qNMR purity in one call, from the raw files to the uncertainty budget."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from copy import deepcopy
from dataclasses import InitVar, dataclass, replace
from math import hypot, isclose, isfinite, sqrt
from numbers import Integral
from os import PathLike, fsdecode
from types import SimpleNamespace

from .acquisition import (
    PILOT_THRESHOLDS, CheckResult, RecoveryEvidence, SignalEvidence, ThresholdProfile, _result, check_recovery,
    check_resolution, check_snr,
)
from .external_calibration import (
    ExternalCalibrationAssay, ExternalCalibrationResult, _calculate_external_calibration, _fill_missing_runs,
    _link_external, _normalize_text, _replace_external_prerequisites,
)
from .integral_checks import _BASELINE_SENSITIVITY_LIMIT, PILOT, IntegralThresholds, _region_thresholds
from .integral_components import ComponentEvidence, ComponentObservation, RegionResult, check_capture_pair, check_region
from .internal_standard import (
    BLANK_SIGNAL_REASON,
    InternalStandardAssay, InternalStandardResult, Prerequisite,
    _calculate_internal_standard, _finite, _has_text, _replace_link_prerequisites, link_spectrum_checks,
)
from .measurement import measure_region
from .processing import (
    QUANTITATIVE_SETTINGS, BaselineShapeAreas, ProcessingSettings, baseline_shape_areas, process,
    uses_quantitative_settings,
)
from .readers import RawFID, UnsupportedInput, read


_ASSESSMENT_TOKEN = object()
_NO_REGION = SimpleNamespace(nominal_area=None)
_FROM_ASSAY = object()


@dataclass(frozen=True)
class AppliedProfiles:
    """The threshold profiles read once at evaluation; the baseline algorithm version is kept apart from its limit."""

    integral: IntegralThresholds
    acquisition: ThresholdProfile
    baseline_limit: float
    baseline_algorithm: str | None
    requested_integral: IntegralThresholds


@dataclass(frozen=True)
class RunEvidence:
    """One run's source and IDs, applied settings, and (role, signal, result or None) regions and checks once measured."""

    run: str
    source: str | None = None
    settings: ProcessingSettings | None = None
    regions: tuple[tuple[str, str, RegionResult | None], ...] | None = None
    acquisition: tuple[CheckResult, ...] | None = None
    run_id: str | None = None
    preparation_id: str | None = None


@dataclass(frozen=True)
class BaselineShapeRegion:
    run: str
    role: str
    signal: str
    bounds_ppm: tuple[float, float] | None
    candidate_areas: BaselineShapeAreas | None


@dataclass(frozen=True)
class BaselineShapeDiagnostic:
    status: str
    reason: str | None = None
    relative_standard_uncertainty_percent: float | None = None
    regions: tuple[BaselineShapeRegion, ...] = ()
    method_version: str = "baseline-shape/1"


def baseline_shape_diagnostic(regions, roles, *, same_tube, quantitative):
    """Summarize fixed baseline candidates for the two sides of an area ratio."""
    if not quantitative:
        return BaselineShapeDiagnostic("not computed", "Requires QUANTITATIVE_SETTINGS")

    def unavailable(reason):
        return BaselineShapeDiagnostic("cannot be assessed", reason, regions=regions)

    grouped = {role: [] for role in roles}
    for region in regions:
        label = f"{region.role} {region.signal}"
        areas = region.candidate_areas
        if areas is None:
            return unavailable(f"{label}: candidate areas unavailable")
        candidates = vars(areas)
        for candidate, area in candidates.items():
            if area is None:
                return unavailable(f"{label}: missing accepted baseline anchor on at least one side")
            if not isfinite(area) or area <= 0:
                return unavailable(f"{label}: {candidate} candidate area must be finite and positive")
        grouped[region.role].append(candidates.values())
    for role, rows in grouped.items():
        if not rows:
            return unavailable(f"{role}: no processed regions for baseline shape")
    first, second = (tuple(map(sum, zip(*rows))) for rows in grouped.values())
    series = (tuple(a / b for a, b in zip(first, second)),) if same_tube else (first, second)
    if any(not isfinite(value) or value <= 0 for values in series for value in values):
        return unavailable("Candidate sums or ratios must be finite and positive")
    terms = (max(abs(value / nominal - 1) for value in alternatives) / sqrt(3)
             for nominal, *alternatives in series)
    percent = 100 * hypot(*terms)
    if not isfinite(percent):
        return unavailable("Relative standard uncertainty is nonfinite")
    return BaselineShapeDiagnostic("computed", relative_standard_uncertainty_percent=percent, regions=regions)


@dataclass(frozen=True)
class DeclaredSignal:
    """A predeclared integral region, its sourced component evidence and its sourced T1."""

    name: str
    bounds_ppm: tuple[float, float] | None
    evidence: ComponentEvidence | None
    t1_s: float | None = None
    t1_source: str | None = None


class _MissingRunId(Prerequisite):
    @property
    def message(self) -> str:
        return f"{self.name} {self.reason}"


@dataclass(frozen=True)
class InternalStandardAssessment:
    """The record linked to the measured evidence and the result calculated from it."""

    _assay: InternalStandardAssay
    _result: InternalStandardResult
    _token: InitVar[object | None] = None
    _profiles: AppliedProfiles | None = None
    _runs: tuple[RunEvidence, ...] = ()
    capture_pairs: tuple[tuple[str, str, ComponentObservation], ...] = ()
    baseline_shape: BaselineShapeDiagnostic | None = None

    def __post_init__(self, _token):
        if _token is not _ASSESSMENT_TOKEN:
            raise TypeError("Use assess_internal_standard to create an assessment")

    @property
    def assay(self) -> InternalStandardAssay:
        return deepcopy(self._assay)

    @property
    def result(self) -> InternalStandardResult:
        return deepcopy(self._result)

    @property
    def profiles(self) -> AppliedProfiles:
        return self._profiles

    @property
    def runs(self) -> tuple[RunEvidence, ...]:
        return deepcopy(self._runs)


@dataclass(frozen=True)
class ExternalCalibrationAssessment:
    """Both runs' record linked to their measured evidence and the result calculated from it."""

    _assay: ExternalCalibrationAssay
    _result: ExternalCalibrationResult
    _token: InitVar[object | None] = None
    _profiles: AppliedProfiles | None = None
    _runs: tuple[RunEvidence, ...] = ()
    capture_pairs: tuple[tuple[str, str, ComponentObservation], ...] = ()
    baseline_shape: BaselineShapeDiagnostic | None = None

    def __post_init__(self, _token):
        if _token is not _ASSESSMENT_TOKEN:
            raise TypeError("Use assess_external_calibration to create an assessment")

    profiles = InternalStandardAssessment.profiles
    runs = InternalStandardAssessment.runs

    @property
    def assay(self) -> ExternalCalibrationAssay:
        return deepcopy(self._assay)

    @property
    def entered_assay(self) -> ExternalCalibrationAssay:
        """Return the assay for editing header observations before reassessment."""
        assay = self.assay
        return replace(assay, _header_values=(), _linked_runs=(),
                       **{role: entered for role, entered, _ in assay._linked_runs})

    @property
    def result(self) -> ExternalCalibrationResult:
        return deepcopy(self._result)


def assess_internal_standard(fid: str | PathLike | RawFID, assay: InternalStandardAssay, *,
                             analyte: Iterable[DeclaredSignal], standard: Iterable[DeclaredSignal],
                             noise_bounds_ppm: tuple[float, float] | None,
                             recovery: RecoveryEvidence,
                             reference: tuple[float, str] | None = None, run_id: str | None = _FROM_ASSAY,
                             preparation_id: str | None = _FROM_ASSAY) -> InternalStandardAssessment:
    """Read, process, measure, check, link and calculate; the only route to a checked pass.

    The assay supplies the quantities and signal names. Measured evidence replaces any entered
    areas, area sources and area unit, and new checks replace earlier generated checks.
    Reference takes a shift in ppm and its source. Use referenced coordinates for all regions.
    A run_id is required for checked purity. If no preparation_id is given, the assay's shared tube id is used.
    The reader's UnsupportedInput exception propagates to the caller, including for a missing file.
    """
    tube = SimpleNamespace(preparation_id=assay.analyte_tube_id if assay.analyte_tube_id == assay.standard_tube_id
                           else None)
    ids = _run_ids(tube, run_id, preparation_id)
    record = _replace_link_prerequisites(assay, (), analyte_areas=None, standard_areas=None,
                                         analyte_area_sources=None, standard_area_sources=None, area_unit=None)
    record, analyte, standard, noise_bounds_ppm, recovery, ids = _normalize_text(deepcopy(
        (record, tuple(analyte), tuple(standard), noise_bounds_ppm, recovery, ids)))
    settings, measurement_checks = _measurement_settings(QUANTITATIVE_SETTINGS, noise_bounds_ppm, reference)
    profiles = _applied_profiles(settings)
    run = RunEvidence("run", **ids)
    id_checks = () if _has_text(run.run_id) else (_MissingRunId("run_id", "cannot be assessed", "is missing"),)
    try:
        raw = _read_raw(fid)
    except UnsupportedInput:
        raise
    except (OSError, ValueError) as error:
        return _withheld(record, profiles, run,
                         Prerequisite("reader", "cannot be assessed", str(error)), *measurement_checks, *id_checks)
    run = replace(run, source=raw.source)
    mismatches = (*_signal_name_checks("analyte", record.analyte_signals, analyte),
                  *_signal_name_checks("standard", record.standard_signals, standard))
    if not mismatches:
        analyte = sorted(analyte, key=lambda signal: record.analyte_signals.index(signal.name))
        standard = sorted(standard, key=lambda signal: record.standard_signals.index(signal.name))
    declared = (*(("analyte", signal) for signal in analyte), *(("standard", signal) for signal in standard))
    declaration, bounds = _declaration_checks(declared)
    declaration = (*id_checks, *declaration)
    input_checks = (*_input_checks(raw), *measurement_checks)
    if mismatches or input_checks:
        return _withheld(record, profiles, run, *declaration, *mismatches, *input_checks)
    measured = _measure(raw, declared, bounds, noise_bounds_ppm, recovery, settings, profiles)
    if isinstance(measured, Prerequisite):
        return _withheld(record, profiles, run, measured, *id_checks)
    regions, failures, acquisition, applied, shape_regions = measured
    run = replace(run, settings=applied, regions=_found(declared, regions), acquisition=tuple(acquisition))
    shape = baseline_shape_diagnostic(shape_regions, ("analyte", "standard"), same_tube=True,
                                      quantitative=uses_quantitative_settings(settings))
    return _link_internal(record, raw, declaration, regions, failures, acquisition, profiles, run, shape)


def _run_ids(run, run_id, preparation_id):
    """Take each ID argument not given from the assay run."""
    return {name: getattr(run, name, None) if value is _FROM_ASSAY else value
            for name, value in (("run_id", run_id), ("preparation_id", preparation_id))}


def _applied_profiles(settings):
    locked = settings.smooth_baseline and settings.lock_spike_anchors
    algorithm = None
    if settings.smooth_baseline or settings.baseline_windows_ppm:
        algorithm = "locked-anchors/1" if locked else "refit-anchors/1"
    requested = PILOT
    return AppliedProfiles(_region_thresholds(requested, locked), PILOT_THRESHOLDS,
                           _BASELINE_SENSITIVITY_LIMIT, algorithm, requested)


def _found(declared, regions):
    remaining = {role: iter(entries) for role, entries in regions.items()}
    found = ((role, signal.name, next(remaining[role])[1]) for role, signal in declared)
    return tuple((role, name, None if region is _NO_REGION else region) for role, name, region in found)


def _measurement_settings(settings, noise_bounds_ppm, reference, label=""):
    checks = []
    if noise_bounds_ppm is not None and (not isinstance(noise_bounds_ppm, tuple) or len(noise_bounds_ppm) != 2
                                         or not all(map(_finite, noise_bounds_ppm))):
        checks.append(Prerequisite(f"{label}noise bounds", "cannot be assessed", "Declare two finite ppm bounds"))
    if reference is not None:
        reference = _normalize_text(reference)
        if (not isinstance(reference, tuple) or len(reference) != 2 or not _finite(reference[0])
                or not _has_text(reference[1])):
            checks.append(Prerequisite(f"{label}reference", "cannot be assessed", "Declare shift_ppm with its source"))
        else:
            settings = replace(settings, reference_shift_ppm=reference[0], reference_source=reference[1])
    return settings, checks


def _measure(raw, declared, bounds, noise_bounds_ppm, recovery, settings, profiles, label=""):
    """Process one run and return its regions, checks and settings, or its processing failure."""
    try:
        spectrum = process(raw, settings, regions_ppm=tuple(span for span in bounds if span),
                           baseline_limit=profiles.baseline_limit)
    except UnsupportedInput:
        raise
    except (ArithmeticError, ValueError) as error:
        return Prerequisite(f"{label}processing", "cannot be assessed", str(error))

    # Every declared signal keeps its region number; a signal without a region result holds a placeholder.
    # A region whose measurement failed keeps its checks but no area.
    regions, signals, failures = {role: [] for role, _ in declared}, [], []
    for (role, signal), span in zip(declared, bounds):
        regions[role].append(("", _NO_REGION))
        if span is None:
            continue
        try:
            region = check_region(spectrum, span, signal.evidence or ComponentEvidence(), raw=raw,
                                  noise_bounds_ppm=noise_bounds_ppm, thresholds=profiles.integral)
            regions[role][-1] = (signal.name, region)
            measured = measure_region(spectrum, signal.name, span, noise_bounds_ppm).to_signal_measurement()
        except ArithmeticError as error:
            failures.append(Prerequisite(f"{signal.name} {role} processing", "cannot be assessed", str(error)))
            name, region = regions[role][-1]
            regions[role][-1] = (name, region if region is _NO_REGION else replace(region, nominal_area=None))
            continue
        signals.append(SignalEvidence(signal.name, signal.t1_s, signal.t1_source, measured))
    extra = recovery.extra_recovery_s
    extra_reason = ("Must be finite and nonnegative" if not _finite(extra) or extra < 0 else
                    "extra_recovery_source is missing"
                    if extra > 0 and not _has_text(recovery.extra_recovery_source) else None)
    evidence = replace(recovery, extra_recovery_s=0.0) if extra_reason else recovery
    # The header check accepted the supplied time; recovery uses the time the sampling grid defines.
    timing = replace(raw, acquisition_time_s=int(raw.complex_points) / float(raw.spectral_width_hz),
                     source_fields={**raw.source_fields, "acquisition_time_s": "complex_points / spectral_width_hz"})
    acquisition, limits = [], profiles.acquisition
    for signal in signals:
        for check, run, limit, unit in (
            ("recovery", lambda: check_recovery(timing, signal, evidence, limits), limits.max_recovery_loss,
             "fraction"),
            ("snr", lambda: check_snr(timing, signal, limits), limits.min_snr, "ratio"),
            ("resolution", lambda: check_resolution(timing, signal, limits), limits.min_points_per_fwhm,
             "points / FWHM"),
        ):
            try:
                acquisition.append(run())
            except ArithmeticError as error:
                acquisition.append(_result(timing, signal, check, limit, unit, (), reason=str(error)))
    if extra_reason:
        failures.append(Prerequisite(f"{label}extra_recovery_s", "cannot be assessed", extra_reason))
        acquisition = [replace(check, status="cannot be assessed", value=None, reason=extra_reason)
                       if check.check == "recovery" and check.value is not None else check for check in acquisition]
    shape_regions = tuple(
        BaselineShapeRegion(label.strip() or "run", role, signal.name, span,
                            baseline_shape_areas(spectrum, span) if span else None)
        for (role, signal), span in zip(declared, bounds)
    )
    return regions, failures, acquisition, spectrum.settings, shape_regions


def _link_internal(record, raw, declaration, regions, failures, acquisition, profiles, run, shape):
    """Link the same-tube run's measured evidence to the record, then calculate."""
    linked = link_spectrum_checks(record, run=str(raw.source), analyte_regions=regions["analyte"],
                                  standard_regions=regions["standard"], acquisition=acquisition)
    # The caller checked supplied names, so a blank region name can only be a placeholder.
    generated = tuple(check for check in linked._linked_prerequisites if check.reason != BLANK_SIGNAL_REASON)
    pairs = _capture_pairs(regions["analyte"], regions["standard"], profiles.integral)
    linked = _replace_link_prerequisites(linked, (*declaration, *failures, *generated, *_pair_prerequisites(pairs)))
    calculation = _calculate_internal_standard(linked)
    status = "cannot be assessed" if calculation.reasons else "pass"
    result = InternalStandardResult(linked.assay_id, status, *calculation)
    return InternalStandardAssessment(linked, result, _ASSESSMENT_TOKEN, profiles, (run,), pairs, shape)


def assess_external_calibration(sample_fid: str | PathLike | RawFID, calibrant_fid: str | PathLike | RawFID,
                                assay: ExternalCalibrationAssay, *, sample: Iterable[DeclaredSignal],
                                calibrant: Iterable[DeclaredSignal],
                                sample_noise_bounds_ppm: tuple[float, float] | None,
                                calibrant_noise_bounds_ppm: tuple[float, float] | None,
                                sample_recovery: RecoveryEvidence,
                                calibrant_recovery: RecoveryEvidence,
                                sample_reference: tuple[float, str] | None = None,
                                calibrant_reference: tuple[float, str] | None = None,
                                sample_run_id: str | None = _FROM_ASSAY,
                                sample_preparation_id: str | None = _FROM_ASSAY,
                                calibrant_run_id: str | None = _FROM_ASSAY,
                                calibrant_preparation_id: str | None = _FROM_ASSAY) -> ExternalCalibrationAssessment:
    """Read, process, measure, check and link each run, then calculate; the only route to a checked pass.

    A failed run blocks purity and leaves the other run's results visible. Measured areas replace entered
    areas; header values and recovery checks come from the current files, never from an earlier link.
    Each run can supply a reference shift in ppm and its source, with regions on that referenced axis.
    Each run's record run_id and preparation_id label it in the report; an ID not given comes from the assay run.
    A missing run record stands as an empty run, so its absence is a named reason.
    UnsupportedInput from either run's reader propagates to the caller, including for a missing file.
    """
    assay, missing_runs = _fill_missing_runs(assay)
    record, sample_inputs, calibrant_inputs, ids = _normalize_text(deepcopy((
        _replace_external_prerequisites(assay, ()),
        (tuple(sample), sample_noise_bounds_ppm, sample_recovery),
        (tuple(calibrant), calibrant_noise_bounds_ppm, calibrant_recovery),
        {"sample": _run_ids(assay.sample, sample_run_id, sample_preparation_id),
         "calibrant": _run_ids(assay.calibrant, calibrant_run_id, calibrant_preparation_id)},
    )))
    base_settings = QUANTITATIVE_SETTINGS
    checks = {role: [] for role in ("sample", "calibrant")}
    evidence, profiles = dict.fromkeys(("sample", "calibrant")), _applied_profiles(base_settings)
    runs = {role: RunEvidence(role, **ids[role]) for role in evidence}
    shape_regions = []
    for role, fid, (signals, noise_bounds_ppm, recovery), reference in (
        ("sample", sample_fid, sample_inputs, sample_reference),
        ("calibrant", calibrant_fid, calibrant_inputs, calibrant_reference),
    ):
        settings, measurement_checks = _measurement_settings(base_settings, noise_bounds_ppm, reference, f"{role} ")
        try:
            raw = _read_raw(fid)
        except UnsupportedInput:
            raise
        except (OSError, ValueError) as error:
            checks[role].append(Prerequisite(f"{role} reader", "cannot be assessed", str(error)))
            checks[role].extend(measurement_checks)
            continue
        evidence[role], runs[role] = (raw, None, ()), replace(runs[role], source=raw.source)
        names = getattr(record, f"{role}_signals")
        mismatches = _signal_name_checks(role, names, signals)
        if not mismatches:
            signals = sorted(signals, key=lambda signal: names.index(signal.name))
        declared = tuple((role, signal) for signal in signals)
        declaration, bounds = _declaration_checks(declared)
        problems = (*mismatches, *(replace(check, name=f"{role} {check.name}") for check in _input_checks(raw)),
                    *measurement_checks)
        checks[role].extend((*declaration, *problems))
        if problems:
            continue
        measured = _measure(raw, declared, bounds, noise_bounds_ppm, recovery, settings, profiles, f"{role} ")
        if isinstance(measured, Prerequisite):
            checks[role].append(measured)
            continue
        regions, failures, acquisition, applied, shape = measured
        shape_regions.extend(shape)
        checks[role].extend(failures)
        evidence[role] = (raw, regions[role], acquisition)
        runs[role] = replace(runs[role], settings=applied, regions=_found(declared, regions),
                             acquisition=tuple(acquisition))
    linked = _link_external(record, evidence, remeasure=True)
    # Supplied names were checked above, so a blank region name can only be a placeholder.
    generated = ((roles, check) for roles, check in zip(linked._linked_check_runs, linked._linked_prerequisites)
                 if check.reason != BLANK_SIGNAL_REASON)
    pairs = _capture_pairs(*(run[1] or () if run else () for run in evidence.values()), profiles.integral)
    linked = _replace_external_prerequisites(linked, (
        *(((role,), check) for role, check in missing_runs.items()),
        *(((role,), check) for role, run_checks in checks.items() for check in run_checks),
        *generated, *((tuple(evidence), check) for check in _pair_prerequisites(pairs)),
    ))
    calculation = _calculate_external_calibration(linked, require_pair_file=True)
    status = "cannot be assessed" if calculation.reasons else "pass"
    result = ExternalCalibrationResult(linked.assay_id, status, *calculation, _assay=linked)
    shape = baseline_shape_diagnostic(tuple(shape_regions), tuple(runs), same_tube=False,
                                      quantitative=uses_quantitative_settings(base_settings))
    return ExternalCalibrationAssessment(linked, result, _ASSESSMENT_TOKEN, profiles, tuple(runs.values()), pairs, shape)


def _capture_pairs(first, second, thresholds):
    """Check the relative target capture of each measured region pair across the ratio."""
    pairs = []
    for first_name, first_region in (entry for entry in first if entry[0]):
        for second_name, second_region in (entry for entry in second if entry[0]):
            pair = check_capture_pair(first_region, second_region, thresholds)
            pairs.append((first_name, second_name, pair))
    return tuple(pairs)


def _pair_prerequisites(pairs):
    return (Prerequisite(f"{first}/{second} capture_pair", pair.status, pair.reason) for first, second, pair in pairs)


def _signal_name_checks(role, expected, supplied):
    status = "cannot be assessed"
    if not isinstance(expected, (tuple, list)) or not expected or any(not _has_text(name) for name in expected):
        return (Prerequisite(f"{role} signals", status, "Declare nonblank signal names"),)
    declared, names = Counter(expected), Counter(signal.name for signal in supplied)
    return (
        *(Prerequisite(f"{name} {role} signal", status, "Signal name is repeated in the assay")
          for name, count in declared.items() if count > 1),
        *(Prerequisite(f"{name} {role} signal", status, "Missing declared signal")
          for name in declared if name not in names),
        *(Prerequisite(f"{name} {role} signal", status, "Undeclared signal") for name in names if name not in declared),
        *(Prerequisite(f"{name} {role} signal", status, f"Supplied {count} times; declared {declared[name]}")
          for name, count in names.items() if name in declared and count != declared[name]),
    )


def _read_raw(fid):
    if fid is None:
        raise ValueError("Missing raw data file")
    raw = fid if isinstance(fid, RawFID) else read(fid)
    try:
        source = fsdecode(raw.source) if isinstance(raw.source, (str, PathLike)) else None
    except UnsupportedInput:
        raise
    except (OSError, TypeError, ValueError):
        source = None
    return _normalize_text(replace(raw, source=source))


def _input_checks(raw):
    status = "cannot be assessed"
    checks = []
    points = raw.complex_points
    if not isinstance(points, Integral) or not _finite(points) or points <= 0:
        checks.append(Prerequisite("complex_points", status, "Must be a positive integer"))
    for name in ("spectral_width_hz", "acquisition_time_s"):
        value = getattr(raw, name)
        if not _finite(value) or value <= 0:
            checks.append(Prerequisite(name, status, "Must be finite and positive"))
    if not checks:
        expected = int(points) / float(raw.spectral_width_hz)
        if not isclose(raw.acquisition_time_s, expected, rel_tol=1e-9):
            reason = (f"{raw.acquisition_time_s:g} s does not match "
                      f"complex_points / spectral_width_hz ({expected:g} s)")
            checks.append(Prerequisite("acquisition_time_s", status, reason))
    if raw.source is None:
        checks.append(Prerequisite("source", status, "Missing raw source"))
    if not _has_text(raw.nucleus) or raw.nucleus.strip().casefold() not in ("1h", "¹h"):
        checks.append(Prerequisite("nucleus", status, f"Expected 1H, got {raw.nucleus!r}"))
    return checks


def _declaration_checks(declared):
    """Block repeated names, unusable bounds and any spectrum interval counted under two names."""
    prerequisites, bounds = [], []
    counts = Counter(signal.name for _, signal in declared)
    for role, signal in declared:
        name = f"{signal.name} {role}"
        if counts[signal.name] > 1:
            prerequisites.append(Prerequisite(f"{name} signal", "cannot be assessed",
                                              "Signal name is declared more than once"))
        span = signal.bounds_ppm
        if not isinstance(span, (tuple, list)) or len(span) != 2 or not all(map(_finite, span)) or span[0] == span[1]:
            prerequisites.append(Prerequisite(f"{name} bounds", "cannot be assessed",
                                              "Declare two distinct finite ppm bounds"))
            span = None
        else:
            span = tuple(sorted(map(float, span)))
            for (other_role, other), other_span in zip(declared, bounds):
                if other_span and span[0] < other_span[1] and other_span[0] < span[1]:
                    prerequisites.append(Prerequisite(f"{name} region", "cannot be assessed",
                                                      f"Overlaps the {other.name} {other_role} region"))
        bounds.append(span)
    return tuple(prerequisites), bounds


def _withheld(record, profiles, run, *checks):
    record = _replace_link_prerequisites(record, checks)
    result = InternalStandardResult(record.assay_id, "cannot be assessed", None, None, None,
                                    reasons=tuple(check.message for check in record.prerequisites if check.status != "pass"))
    shape = baseline_shape_diagnostic((), ("analyte", "standard"), same_tube=True,
                                      quantitative=uses_quantitative_settings(QUANTITATIVE_SETTINGS))
    return InternalStandardAssessment(record, result, _ASSESSMENT_TOKEN, profiles, (run,), baseline_shape=shape)
