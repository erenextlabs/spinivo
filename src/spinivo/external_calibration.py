"""Calculate external-calibration qNMR purity and its GUM uncertainty budget."""

from __future__ import annotations

import re
from collections.abc import Iterable
from copy import deepcopy
from dataclasses import dataclass, field, fields, is_dataclass, replace
from fractions import Fraction
from math import isclose
from numbers import Rational
from os import stat, walk
from os.path import abspath, expanduser, join, samefile
from pathlib import Path
from stat import S_ISREG
from typing import NamedTuple, TYPE_CHECKING

from .internal_standard import (
    CALCULATION_ONLY, Correlation, Covariance, CoverageFactor, Prerequisite, Quantity, UncertaintyBudget,
    _FRACTION_UNITS, _area_sum, _finite, _link_checks, _link_regions, _link_units, _ratio_quantity,
    _replace_link_prerequisites, _uncertainty_budget, _validate_correlations,
    _validate_number as _validate_finite_number,
)

if TYPE_CHECKING:
    from .acquisition import CheckResult
    from .integral_components import RegionResult
    from .readers.model import RawFID


@dataclass(frozen=True)
class ExternalRun:
    """Record area accumulation and nominal settings, with temperature in Celsius.

    reader_scale is the multiplier the reader applied to the stored intensity (1 for stored values).
    Supply reader_scale_source for independent evidence applicable to this run when linking.
    D1, AQ, SW, points and frequency are recorded per run; only the frequency is compared, within tolerance.
    """

    run_id: str | None
    preparation_id: str | None
    area: float | None
    scans: int | None
    area_mode: str | None
    receiver_gain: float | None
    probe: str | None
    pulse_sequence: str | None
    decoupling: str | None
    temperature: float | None
    solvent: str | None
    tube: str | None
    instrument: str | None = None
    excitation: str | None = None
    fill: str | None = None
    reader_scale: float | None = None
    relaxation_delay_s: float | None = None
    acquisition_time_s: float | None = None
    spectral_width_hz: float | None = None
    complex_points: int | None = None
    observe_frequency_mhz: float | None = None
    storage_scale: str | None = None
    calibrated_flip_angle_deg: float | None = None
    reader_scale_source: str | None = None


@dataclass(frozen=True)
class PairRecord:
    """Tie both runs and preparations to a record or analyst attestation.

    Relative run and source paths share the working directory as their base.
    """

    kind: str | None
    sample_run_id: str | None
    sample_preparation_id: str | None
    calibrant_run_id: str | None
    calibrant_preparation_id: str | None
    source: str | None


@dataclass(frozen=True)
class ExternalCalibrationAssay:
    """K is calibrant response divided by sample response. D is drift. Both stay separate from S.

    S = (sample stored-to-common multiplier) / (calibrant multiplier) for the per-scan stored area ratio.
    S_gains and S_storage_scales identify the settings S.source covers in sample, calibrant order.
    K_covers lists condition names whose differences are justified by K.source. TRANSFER_TOLERANCES gives
    temperature in Celsius and observe frequency as a relative difference.
    Declare c_s_route as "certified" (c_s in mmol/L) or "prepared" (four preparation quantities).
    Declare the signal names for each region when linking measured areas.
    """

    assay_id: str | None
    sample: ExternalRun | None
    calibrant: ExternalRun | None
    N_a: int | None
    N_s: int | None
    M_a: Quantity | None
    V_a: Quantity | None
    m_a: Quantity | None
    K: Quantity | None
    D: Quantity | None
    R_relative_u: float | None
    R_u_source: str
    c_s: Quantity | None = None
    m_s: Quantity | None = None
    w_s: Quantity | None = None
    M_s: Quantity | None = None
    V_s: Quantity | None = None
    K_covers: tuple[str, ...] = ()
    correlations: tuple[Correlation | Covariance, ...] = ()
    k: CoverageFactor | None = None
    prerequisites: tuple[Prerequisite, ...] = ()
    c_s_route: str | None = None
    R_source_ids: tuple[str, ...] = ()
    sample_signals: tuple[str, ...] = ()
    calibrant_signals: tuple[str, ...] = ()
    pair: PairRecord | None = None
    S: Quantity | None = None
    S_gains: tuple[float, float] | None = None
    S_storage_scales: tuple[str, str] | None = None
    _linked_prerequisites: tuple[Prerequisite, ...] = field(default=(), repr=False, compare=False)
    # Run roles each generated prerequisite depends on, in prerequisite order.
    _linked_check_runs: tuple[tuple[str, ...], ...] = field(default=(), repr=False, compare=False)
    _header_values: tuple[tuple[str, str, object], ...] = field(default=(), repr=False, compare=False)
    _linked_runs: tuple[tuple[str, ExternalRun, ExternalRun], ...] = field(default=(), repr=False, compare=False)


@dataclass(frozen=True)
class ExternalCalibrationResult:
    """Concentrations use mmol/L for mg, mL and g/mol, or the paired units as supplied."""

    assay_id: str | None
    status: str
    J_a: float | None
    J_s: float | None
    R: float | None
    value: float | None = None
    budget: UncertaintyBudget | None = None
    reasons: tuple[str, ...] = ()
    c_s: float | None = None
    c_a: float | None = None
    concentration_unit: str | None = None
    transfer_tolerances: tuple[tuple[str, float], ...] = ()
    _assay: ExternalCalibrationAssay | None = field(default=None, repr=False, compare=False)


class _ExternalCalibrationCalculation(NamedTuple):
    J_a: float | None
    J_s: float | None
    R: float | None
    value: float | None
    budget: UncertaintyBudget | None
    reasons: tuple[str, ...]
    c_s: float | None
    c_a: float | None
    concentration_unit: str | None
    transfer_tolerances: tuple[tuple[str, float], ...]


_PREPARATION = ("m_s", "w_s", "M_s", "V_s")
_UNIT_PAIRS = (("m_s", "m_a", "mg"), ("V_s", "V_a", "mL"), ("M_s", "M_a", "g/mol"))
_CONDITIONS = ("instrument", "probe", "pulse_sequence", "decoupling", "excitation", "temperature", "solvent", "tube",
               "fill", "calibrated_flip_angle_deg")
_RECORDED = ("relaxation_delay_s", "acquisition_time_s", "spectral_width_hz", "complex_points", "observe_frequency_mhz")
TRANSFER_TOLERANCES = {"temperature": 0.1, "observe_frequency_mhz": 1e-4}
_MISSING_TEXT = frozenset(("", "unknown", "not known", "not found", "not available", "not recorded", "not documented",
                           "not supplied", "not measured", "not specified", "not applicable", "missing", "unspecified",
                           "unavailable", "na", "n a", "none", "null", "nan", "tbd", "undetermined", "nf", "n f",
                           "not given", "not provided", "not reported", "not stated", "no data", "no record", "tba"))


def _missing(value):
    """The one missing-value test for evidence text: blank, NF (NOT FOUND), unknown, n/a and similar."""
    return value is None or (isinstance(value, str)
                             and re.sub(r"[\W_]+", " ", value.casefold()).strip() in _MISSING_TEXT)


def _has_text(value):
    return isinstance(value, str) and not _missing(value)


def _normalize_text(value):
    """Replace missing-text markers throughout an input record or collection of records."""
    if isinstance(value, str):
        return None if _missing(value) else value
    if isinstance(value, (tuple, list)):
        return tuple(map(_normalize_text, value))
    if isinstance(value, dict):
        return {key: _normalize_text(item) for key, item in value.items()}
    if is_dataclass(value) and not isinstance(value, type):
        return replace(value, **{f.name: _normalize_text(getattr(value, f.name)) for f in fields(value) if f.init})
    return value


def _validate_number(name, value, reasons, *, zero_ok=False):
    return _validate_finite_number(name, None if _missing(value) else value, reasons, zero_ok=zero_ok)


def link_external_checks(assay: ExternalCalibrationAssay, *, sample_raw: RawFID, calibrant_raw: RawFID,
                         sample_regions: Iterable[tuple[str, RegionResult]],
                         calibrant_regions: Iterable[tuple[str, RegionResult]],
                         sample_acquisition: Iterable[CheckResult],
                         calibrant_acquisition: Iterable[CheckResult]) -> ExternalCalibrationAssay:
    """Attach both runs' nominal areas, header settings and qualification checks."""
    assay, missing_runs = _fill_missing_runs(assay)
    assay = _normalize_text(_replace_external_prerequisites(assay, ()))
    evidence = {role: _normalize_text((raw, tuple(regions), tuple(acquisition)))
                for role, raw, regions, acquisition in (
                    ("sample", sample_raw, sample_regions, sample_acquisition),
                    ("calibrant", calibrant_raw, calibrant_regions, calibrant_acquisition),
                ) if role not in missing_runs}
    return _link_external(assay, evidence, missing_runs=missing_runs)


def _replace_external_prerequisites(assay, scoped_checks, **changes):
    scoped_checks = tuple(scoped_checks)
    return _replace_link_prerequisites(assay, (check for _, check in scoped_checks),
                                       _linked_check_runs=tuple(runs for runs, _ in scoped_checks), **changes)


def _link_external(assay, evidence, *, remeasure=False, missing_runs=None):
    """Record the current headers and link any measured regions, replacing previous header-derived values."""
    runs = {role: getattr(assay, role) for role in ("sample", "calibrant")}
    for role, entered, linked in assay._linked_runs:
        run = runs[role]
        changes = {name: value for name, value in vars(entered).items()
                   if getattr(run, name) == getattr(linked, name)
                   and not (name == "reader_scale" and _has_text(run.reader_scale_source))}
        runs[role] = replace(run, **changes)
    entered_runs = runs.copy()
    if remeasure:
        runs = {role: replace(run, area=None, storage_scale=None, calibrated_flip_angle_deg=None)
                for role, run in runs.items()}
    missing_runs = missing_runs or {}
    scoped = [((role,), check) for role, check in missing_runs.items()]
    region_units, unit_roles, headers = {}, [], []
    for role, run_evidence in evidence.items():
        if run_evidence is not None:
            raw, regions, acquisition = run_evidence
            header = _run_header(raw)
            headers.extend((role, name, value) for name, value in header.items())
            checks, units = [], {}
            runs[role] = _link_run(role, runs[role], raw, regions, getattr(assay, f"{role}_signals"), acquisition,
                                   checks, units, header)
            scoped.extend(((role,), check) for check in checks)
            if units:
                unit_roles.append(role)
                region_units.update(units)
    checks = []
    _link_units(region_units, checks)
    scoped.extend((tuple(unit_roles), check) for check in checks)
    linked_runs = tuple((role, entered_runs[role], run) for role, run in runs.items() if role not in missing_runs)
    # Keep absent records absent on direct links so a retry sees the same missing input.
    runs.update(dict.fromkeys(missing_runs))
    return _replace_external_prerequisites(assay, scoped, _header_values=tuple(headers),
                                           _linked_runs=linked_runs, **runs)


def _link_run(role, run, raw, regions, expected, acquisition, prerequisites, region_units, header):
    source = str(raw.source)
    independent_scale = _has_text(run.reader_scale_source) and not _missing(run.reader_scale)
    run = replace(run, reader_scale=run.reader_scale if independent_scale else None,
                  reader_scale_source=run.reader_scale_source if independent_scale else None,
                  storage_scale=None, calibrated_flip_angle_deg=None)
    if run.run_id != source:
        if not _has_text(run.run_id):
            return run
        prerequisites.append(Prerequisite(f"{role} run_id", "cannot be assessed",
                                          f"record {run.run_id!r} differs from raw source {source!r}"))
        return run
    run = replace(_link_header(role, run, header, prerequisites), storage_scale=header["storage_scale"])
    if regions is None:
        return run
    _, signals = _link_regions(role, regions, expected, source, prerequisites, region_units)
    _link_checks(f"{role} ", signals, source, acquisition, prerequisites)
    angles = [next((item.value for item in check.inputs if item.name == "calibrated_flip_angle_deg"), None)
              for check in acquisition if check.run == source and check.signal in signals and check.check == "recovery"]
    angle = angles[0] if angles and all(_finite(value) and value == angles[0] for value in angles) else None

    area_reasons = []
    area = _linked_area_sum(role, [region.nominal_area for _, region in regions], area_reasons)
    if area is None:
        prerequisites.extend(Prerequisite(f"{role} linked area", "cannot be assessed", reason)
                             for reason in area_reasons)
    return replace(run, area=area, calibrated_flip_angle_deg=angle)


def _linked_area_sum(role, areas, reasons):
    if any(isinstance(area, Fraction) for area in areas) and all(_finite(area) for area in areas):
        # Keep the exact sum through scan normalization.
        total = sum(Fraction(int(area.numerator), int(area.denominator)) if isinstance(area, Rational)
                    else Fraction(float(area)) for area in areas)
        if _finite(total):
            return total
        reasons.append(f"{role} area sum is outside the numeric range")
        return None
    return _area_sum(role, areas, reasons)


def _run_header(raw):
    temperature = round(raw.temperature_k - 273.15, 6) if raw.temperature_k is not None else None
    return {"scans": raw.scans, "receiver_gain": raw.receiver_gain, "reader_scale": raw.intensity_scale,
            "pulse_sequence": raw.pulse_program, "temperature": temperature, "storage_scale": raw.storage_scale,
            **{name: getattr(raw, name) for name in _RECORDED}}


def _link_header(role, run, header, prerequisites):
    updates = {}
    for name, value in header.items():
        if name == "storage_scale":
            continue
        record = getattr(run, name)
        if _missing(value):
            if name != "reader_scale" and not _missing(record):
                prerequisites.append(Prerequisite(f"{role} {name}", "cannot be assessed",
                                                  f"header has no value; entered {record} is not header evidence"))
        elif _missing(record):
            updates[name] = value
        elif record != value:
            prerequisites.append(Prerequisite(f"{role} {name}", "cannot be assessed",
                                              f"record {record} differs from header {value}"))
    return replace(run, **updates)


def _scan_area(run):
    if not _finite(run.area):
        return None
    area = run.area if isinstance(run.area, Rational) else float(run.area)
    if run.area_mode == "averaged":
        return area
    if run.area_mode == "summed" and _finite(run.scans) and run.scans > 0 and not run.scans % 1:
        return area / int(run.scans)
    return None


def _per_scan_area(role, run, reasons):
    for name in ("run_id", "preparation_id"):
        if not _has_text(getattr(run, name)):
            reasons.append(f"{role} {name} is missing")
    scans_ok = _validate_number(f"{role} scans", run.scans, reasons)
    if scans_ok and run.scans % 1:
        reasons.append(f"{role} scans must be a whole count")
    mode = run.area_mode
    if _missing(mode):
        reasons.append(f"{role} area_mode is missing")
    elif mode not in ("summed", "averaged"):
        reasons.append(f"{role} area_mode must be summed or averaged")
    scale_ok = _validate_number(f"{role} reader_scale", run.reader_scale, reasons)
    _validate_number(f"{role} area", run.area, reasons)
    if scale_ok and float(run.reader_scale) == 0:
        reasons.append(f"{role} reader_scale is outside the numeric range")
        return None
    area = _scan_area(run)
    if not scale_ok or area is None:
        return None
    try:
        area = float(area / Fraction(run.reader_scale) if isinstance(area, Rational) else area / run.reader_scale)
    except OverflowError:
        area = None
    if not _finite(area) or run.area != 0 and area == 0:
        reasons.append(f"{role} per-scan area is outside the numeric range")
        return None
    return area


def _compatible(name, a, s, tolerances):
    tolerance = tolerances.get(name)
    if tolerance is None:
        return a == s
    limit = tolerance * (max(abs(a), abs(s)) if name == "observe_frequency_mhz" else 1)
    difference = abs(a - s)
    return difference <= limit or isclose(difference, limit, rel_tol=1e-12)


def _validate_transfer(assay, reasons, tolerances):
    sample, calibrant = assay.sample, assay.calibrant
    if _has_text(sample.run_id) and sample.run_id == calibrant.run_id:
        reasons.append("sample and calibrant run_id must differ")
    covered = assay.K_covers
    if not isinstance(covered, (tuple, list)):
        reasons.append("K_covers must be a tuple or list of condition names")
        covered = ()
    for name in covered:
        if name in ("instrument", "receiver_gain"):
            reasons.append(f"K cannot cover {name}")
        elif name not in _CONDITIONS:
            reasons.append(f"K_covers contains unknown condition {name!r}")
    for name in (*_CONDITIONS, *_RECORDED):
        a, s = getattr(sample, name), getattr(calibrant, name)
        before = len(reasons)
        for role, value in (("sample", a), ("calibrant", s)):
            label = f"{role} {name}"
            if name in _RECORDED or name == "calibrated_flip_angle_deg":
                if _validate_number(label, value, reasons):
                    if name == "complex_points" and value % 1:
                        reasons.append(f"{label} must be a whole count")
                    elif name == "calibrated_flip_angle_deg" and value > 180:
                        reasons.append(f"{label} must be at most 180")
            elif not (_finite(value) if name == "temperature" else _has_text(value)):
                reasons.append(f"{label} {'is missing' if _missing(value) else 'is invalid'}")
        if ((name in _RECORDED and name != "observe_frequency_mhz") or len(reasons) != before
                or _compatible(name, a, s, tolerances)):
            continue
        if name in ("instrument", "observe_frequency_mhz"):
            reasons.append(f"{name} differs; transfer between instruments is not supported")
        elif name not in covered:
            reasons.append(f"{name} differs without K evidence")


def _same_file(path, other):
    if path == other:
        return True
    try:
        return samefile(path, other)
    except OSError:
        return False


def _spectrum_file(path):
    return path.suffix.casefold() == ".jdf" or path.name.casefold() in ("fid", "ser")


def _raise_walk_error(error):
    raise error


def _run_file(source_path, source, path, run):
    """Check both path containment and file identity against the run."""
    if source in run.parents or _same_file(source, run):
        return True
    if path.name.casefold() in ("fid", "ser") and run.is_file():
        path, run = path.parent, run.parent
    if not run.is_dir() and (run.is_file() or _spectrum_file(path) or _spectrum_file(run)):
        return False
    if source.is_relative_to(run) or Path(abspath(source_path)).is_relative_to(abspath(path)):
        return True
    return run.is_dir() and source.is_file() and any(
        _same_file(source, Path(folder, name)) for folder, _, names in walk(run, onerror=_raise_walk_error)
        for name in names
    )


def _check_pair_file(source, reasons):
    try:
        if not S_ISREG(stat(source).st_mode):
            reasons.append("pair source must be a readable file")
            return
        with open(source, "rb") as stream:
            if not stream.read(1):
                reasons.append("pair source file is empty")
    except FileNotFoundError:
        reasons.append("pair source file is missing")
    except OSError as error:
        reasons.append(f"pair source file could not be read: {error}")


def _validate_pair(assay, reasons, *, require_file=False):
    pair = assay.pair
    if not isinstance(pair, PairRecord):
        reasons.append("pair record is missing" if pair is None else "pair must be a PairRecord")
        return
    if pair.kind not in ("record", "attestation"):
        reasons.append("pair kind must be record or attestation")
    if not _has_text(pair.source) and not (isinstance(pair.source, str) and pair.source.strip() in (".", "..")):
        reasons.append("pair source is missing")
    else:
        try:
            base = Path.cwd()
            source_path = join(base, expanduser(pair.source if require_file else pair.source.strip()))
            entered_path = Path(source_path)
            source = entered_path.resolve()
            paths = [base / Path(run.run_id).expanduser() for run in (assay.sample, assay.calibrant)
                     if _has_text(run.run_id)]
            runs = [path.resolve() for path in paths]
            if (source.is_dir() or _spectrum_file(entered_path) or _spectrum_file(source)
                    or any(_run_file(source_path, source, path, run) for path, run in zip(paths, runs))):
                reasons.append("pair source names a run file or folder, not a pair record")
            elif require_file:
                _check_pair_file(source_path, reasons)
            if len(runs) == 2 and assay.sample.run_id != assay.calibrant.run_id and _same_file(*runs):
                reasons.append("sample and calibrant run_id must differ")
        except OSError as error:
            reasons.append(f"pair source or run folder could not be inspected: {error}")
        except (RuntimeError, ValueError):
            reasons.append("pair source or run path could not be resolved")
    for role in ("sample", "calibrant"):
        for name in ("run_id", "preparation_id"):
            value, record = getattr(pair, f"{role}_{name}"), getattr(getattr(assay, role), name)
            if not _has_text(value):
                reasons.append(f"pair {role}_{name} is missing")
            elif value != record:
                reasons.append(f"pair {role}_{name} {value!r} differs from the {role} record {record!r}")


def _validate_scale(assay, reasons):
    gains = (assay.sample.receiver_gain, assay.calibrant.receiver_gain)
    scales = (assay.sample.storage_scale, assay.calibrant.storage_scale)
    for role, gain, scale in zip(("sample", "calibrant"), gains, scales):
        _validate_number(f"{role} receiver_gain", gain, reasons)
        if _missing(scale):
            reasons.append(f"{role} storage_scale is missing")
        elif not isinstance(scale, str) or not (scale == "JEOL" or re.fullmatch(r"NC=-?\d+", scale)):
            reasons.append(f"{role} storage_scale is invalid")
    if _missing(assay.S):
        reasons.append("equal nominal receiver gains are not scale evidence" if gains[0] == gains[1]
                       else "receiver gains differ without a sourced S correction")
        return
    stated = assay.S_gains
    if _missing(stated):
        reasons.append("S_gains is missing")
    elif not isinstance(stated, (tuple, list)) or not all(map(_finite, stated)) or tuple(stated) != gains:
        reasons.append(f"S_gains {stated!r} do not cover the run receiver gains {gains!r}")
    stated = assay.S_storage_scales
    if _missing(stated):
        reasons.append("S_storage_scales is missing")
    elif not isinstance(stated, (tuple, list)) or tuple(stated) != scales:
        reasons.append(f"S_storage_scales {stated!r} do not cover the run storage scales {scales!r}")


def _calibrant_concentration(assay, reasons):
    if assay.c_s_route == "certified":
        return assay.c_s.value, "mmol/L"
    m, w, M, V = (getattr(assay, name) for name in _PREPARATION)
    use_mmol = (m.unit, V.unit, M.unit) == ("mg", "mL", "g/mol")
    unit = "mmol/L" if use_mmol else f"{m.unit} / (({M.unit}) * {V.unit})"
    try:
        c_s = (m.value * w.value) / (M.value * V.value) * (1000 if use_mmol else 1)
        if _finite(c_s) and c_s > 0:
            return c_s, unit
    except ArithmeticError:
        pass
    reasons.append("calculated calibrant concentration is outside the numeric range")
    return None, None


def _fill_missing_runs(assay):
    """Supply empty fields for absent runs and discard generated evidence that used them."""
    runs = {role: ExternalRun(**dict.fromkeys(field.name for field in fields(ExternalRun) if field.init))
            for role in ("sample", "calibrant") if getattr(assay, role) is None}
    if not runs:
        return assay, {}
    checks = {role: Prerequisite(f"{role} run", "cannot be assessed", "Missing run record") for role in runs}
    scoped = ((roles, check) for roles, check in zip(assay._linked_check_runs, assay._linked_prerequisites)
              if runs.keys().isdisjoint(roles))
    return _replace_external_prerequisites(
        assay, scoped, _header_values=tuple(row for row in assay._header_values if row[0] not in runs),
        _linked_runs=tuple(row for row in assay._linked_runs if row[0] not in runs), **runs,
    ), checks


def calculate_external_calibration(assay: ExternalCalibrationAssay) -> ExternalCalibrationResult:
    """Preserve computable per-scan areas and their ratio when purity is blocked.

    A direct calculation never passes; purity is reported as calculation only.
    Only spinivo.pipeline.assess_external_calibration can pass.
    """
    assay, missing_runs = _fill_missing_runs(deepcopy(assay))
    if missing_runs:
        assay = replace(assay, prerequisites=(*assay.prerequisites, *missing_runs.values()))
    calculation = _calculate_external_calibration(replace(assay, correlations=_normalize_text(assay.correlations)))
    calculation = calculation._replace(reasons=calculation.reasons or (CALCULATION_ONLY,))
    return ExternalCalibrationResult(assay.assay_id, "cannot be assessed", *calculation, _assay=assay)


def _calculate_external_calibration(assay, *, require_pair_file=False) -> _ExternalCalibrationCalculation:
    """Return computable values and reasons blocking purity."""
    tolerances = TRANSFER_TOLERANCES.copy()
    applied = tuple(tolerances.items())
    reasons = []
    J_a = _per_scan_area("sample", assay.sample, reasons)
    J_s = _per_scan_area("calibrant", assay.calibrant, reasons)
    R = J_a / J_s if J_a is not None and J_s else None
    if R is not None:
        if not _finite(R):
            reasons.append("R is outside the numeric range")
            R = None
        elif R == 0 and J_a != 0:
            reasons.append("R is too small to represent")
    _validate_transfer(assay, reasons, tolerances)
    _validate_pair(assay, reasons, require_file=require_pair_file)
    _validate_scale(assay, reasons)
    for name in ("assay_id", "R_u_source"):
        if not _has_text(getattr(assay, name)):
            reasons.append(f"{name} is missing")
    for name in ("N_a", "N_s"):
        count = getattr(assay, name)
        if _validate_number(name, count, reasons) and count % 1:
            reasons.append(f"{name} must be a whole count")
    _validate_number("R_relative_u", assay.R_relative_u, reasons, zero_ok=True)

    route = assay.c_s_route
    names = ("c_s",) if route == "certified" else _PREPARATION if route == "prepared" else ()
    inactive = _PREPARATION if route == "certified" else ("c_s",)
    collision = any(getattr(assay, name) is not None for name in inactive)
    if not names:
        reasons.append("c_s_route must be declared as certified or prepared")
    elif collision:
        reasons.append(f"only the declared {route} c_s inputs may be supplied")
    quantities = {name: getattr(assay, name) for name in (*names, "M_a", "V_a", "m_a", "K", "D", "S")}
    invalid = set()
    for name, quantity in quantities.items():
        before = len(reasons)
        if not isinstance(quantity, Quantity):
            reasons.append(f"{name} is missing" if _missing(quantity) else f"{name} must be a Quantity")
            invalid.add(name)
            quantities[name] = None
            continue
        if _validate_number(f"{name}.value", quantity.value, reasons) and name == "w_s" and quantity.value > 1:
            reasons.append("w_s.value must be at most 1")
        if not _has_text(quantity.unit):
            reasons.append(f"{name}.unit is missing")
        elif name == "w_s" and quantity.unit not in _FRACTION_UNITS:
            reasons.append(f"w_s.unit must be one of {', '.join(_FRACTION_UNITS)}")
        elif name in ("K", "D", "S") and quantity.unit not in ("1", "ratio", "fraction"):
            reasons.append(f"{name}.unit must be dimensionless")
        elif name == "c_s" and quantity.unit != "mmol/L":
            reasons.append("c_s.unit must be mmol/L")
        if not _has_text(quantity.source):
            reasons.append(f"{name}.source is missing")
        if len(reasons) != before:
            invalid.add(name)
        _validate_number(f"{name}.standard_u", quantity.standard_u, reasons, zero_ok=name != "S")
    c_s, concentration_unit = (None, None)
    if names and not collision and not invalid.intersection(names):
        c_s, concentration_unit = _calibrant_concentration(assay, reasons)
    for standard, sample, unit in _UNIT_PAIRS:
        a, s = quantities[sample], quantities.get(standard)
        if route == "certified" and a is not None and a.unit != unit:
            reasons.append(f"{sample}.unit must be {unit} for certified c_s")
        if s is not None and a is not None and _has_text(s.unit) and _has_text(a.unit) and s.unit != a.unit:
            reasons.append(f"{standard} and {sample} must use the same unit")

    quantities = {"R": _ratio_quantity(R, assay), **quantities}
    correlations = _validate_correlations(assay.correlations, quantities, reasons)
    if assay.k is not None:
        _validate_number("k.value", assay.k.value, reasons)
        if not _has_text(assay.k.source):
            reasons.append("k.source is missing")
    reasons.extend(f"{check.name}: {check.reason}" for check in assay.prerequisites if check.status != "pass")
    if not reasons:
        try:
            m_a, V_a, M_a = assay.m_a, assay.V_a, assay.M_a
            scale = 1000 if concentration_unit == "mmol/L" else 1
            before_K = R * assay.S.value * (assay.N_s / assay.N_a) * c_s * assay.D.value
            sample_factor = (M_a.value * V_a.value / m_a.value) / scale
            c_a = before_K * assay.K.value
            value = (before_K * sample_factor) * assay.K.value
            if not all(_finite(item) and item > 0 for item in (c_a, value)):
                raise ValueError("calculated concentration or purity is outside the numeric range")
            budget = _uncertainty_budget(value, quantities, correlations, assay.k, ("M_s", "V_s", "m_a"))
        except (ArithmeticError, ValueError) as error:
            reasons.append(f"external calibration calculation: {error}")
        else:
            return _ExternalCalibrationCalculation(J_a, J_s, R, value, budget, (), c_s, c_a, concentration_unit, applied)
    return _ExternalCalibrationCalculation(J_a, J_s, R, None, None, tuple(reasons), c_s, None, concentration_unit, applied)
