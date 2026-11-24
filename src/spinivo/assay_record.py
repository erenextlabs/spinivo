"""Load one versioned JSON assay record (spinivo.assay/1) and assess it through the route its mode names."""

from __future__ import annotations

import json
from dataclasses import MISSING, dataclass, fields, replace
from os import PathLike
from os.path import join
from pathlib import Path

from .acquisition import RecoveryEvidence
from .external_calibration import _RECORDED, ExternalCalibrationAssay, ExternalRun, PairRecord, _missing
from .integral_components import ComponentEvidence, Line
from .internal_standard import (
    Correlation, Covariance, CoverageFactor, InternalStandardAssay, Prerequisite, Quantity, _finite, _has_text,
    calculate_internal_standard,
)
from .pipeline import DeclaredSignal, assess_external_calibration, assess_internal_standard

SCHEMA = "spinivo.assay/1"
MODES = ("internal_standard", "external_calibration", "calculation_only")
_ROLES = {"internal_standard": ("analyte", "standard"), "external_calibration": ("sample", "calibrant")}
_QUANTITIES = {"internal_standard": ("M_a", "M_s", "m_a", "m_s", "w_s"),
               "external_calibration": ("M_a", "V_a", "m_a", "K", "D", "c_s", "m_s", "w_s", "M_s", "V_s")}


class RecordError(ValueError):
    """The record is not valid JSON or breaks the spinivo.assay/1 layout; the message starts with the field path."""


@dataclass(frozen=True)
class AssayRecord:
    """The assay, each run's checked entries keyed "run" or by role, and the pipeline keyword inputs, if any."""

    mode: str
    assay: InternalStandardAssay | ExternalCalibrationAssay
    runs: dict[str, dict]
    inputs: dict[str, object]


def _fail(path, problem):
    raise RecordError(f"{path or 'record'}: {problem}")


def _scalar(kinds, label):
    def convert(value, path):
        if value is None or isinstance(value, kinds) and not isinstance(value, bool):
            return value
        _fail(path, f"must be {label} or null")
    return convert


_numeric, _text = _scalar((int, float), "a number"), _scalar(str, "text")


def _number(value, path):
    value = _numeric(value, path)
    if value is not None and not _finite(value):
        _fail(path, "must be finite and within the float range")
    return value


def _json_float(token):
    value = float(token)
    # Keep nonzero values that underflow to zero distinguishable during field validation.
    if value == 0 and token.lower().split("e", 1)[0].strip("0.-"):
        return float("nan")
    return value


def _json_int(token):
    try:
        return int(token)
    except ValueError:
        return float(token)


class _Members(dict):
    """A decoded JSON object that keeps the first repeated member name."""

    repeated = None


def _json_object(pairs):
    members = _Members()
    for key, value in pairs:
        if members.repeated is None and key in members:
            members.repeated = key
        members[key] = value
    return members


def _list(item, size=None):
    def convert(value, path):
        if value is None:
            return None if size else ()
        if not isinstance(value, list) or size and len(value) != size:
            _fail(path, f"must be a list of {size} items or null" if size else "must be a list or null")
        return tuple(item(entry, f"{path}[{index}]") for index, entry in enumerate(value))
    return convert


def _object(spec, cls=None, *, allow_null=True, **absent):
    """Convert an object's entries and reject unknown keys; build cls, with absent fields null unless given."""
    def convert(value, path):
        if value is None and allow_null:
            return None if cls else {}
        if not isinstance(value, dict):
            _fail(path, "must be an object or null" if allow_null else "must be an object")
        prefix = f"{path}." if path else ""
        key = getattr(value, "repeated", None)
        if key is not None:
            _fail(f"{prefix}{key}", "repeated field")
        for key in value:
            if key not in spec:
                _fail(f"{prefix}{key}", "unknown field")
        values = {key: spec[key](entry, f"{prefix}{key}") for key, entry in value.items()}
        return _build(cls, {**absent, **values}, path) if cls else values
    return convert


def _build(cls, values, path):
    nulls = {f.name: None for f in fields(cls) if f.default is MISSING and f.default_factory is MISSING}
    try:
        return cls(**{**nulls, **values})
    except (TypeError, ValueError) as error:
        _fail(path, str(error))


def _texts(*names):
    return dict.fromkeys(names, _text)


_QUANTITY = _object({"value": _number, "standard_u": _number, **_texts("unit", "source"), "source_ids": _list(_text)},
                    Quantity)
_LINE_FIELDS = ("ppm", "area", "fwhm_hz")
_LINE = _object(dict.fromkeys(_LINE_FIELDS, _number), allow_null=False)
_EVIDENCE = _object({"target_lines": _list(_LINE), "foreign_lines": _list(_LINE), "coupling_hz": _number,
                     "satellite_fraction": _number, **_texts("target_source", "foreign_source", "satellite_convention",
                                                             "satellite_source")})
_RECOVERY = _object({**dict.fromkeys(("calibrated_flip_angle_deg", "calibrated_90_us", "extra_recovery_s"), _number),
                     **_texts("calibration_source", "extra_recovery_source")}, RecoveryEvidence, extra_recovery_s=None)
_INDEPENDENCE = {"source": _text, "independence_reason": _text}
_PASSED = ("assay_id", "N_a", "N_s", "R_relative_u", "R_u_source", "R_source_ids", "k", "c_s_route", "K_covers")
_SPEC = {
    **_texts("schema", "mode", "notes", "assay_id", "R_u_source"), "R_source_ids": _list(_text),
    **dict.fromkeys(("N_a", "N_s", "R_relative_u"), _number),
    "correlations": _list(_object({"inputs": _list(_text, 2), "r": _number, **_INDEPENDENCE}, Correlation)),
    "covariances": _list(_object({"inputs": _list(_text, 2), "value": _number, "units": _list(_text, 2),
                                  **_INDEPENDENCE}, Covariance)),
    "k": _object({"value": _number, "source": _text}, CoverageFactor),
}
_EXTERNAL = {
    "c_s_route": _text, "K_covers": _list(_text),
    "pair": _object(_texts("kind", "sample_run_id", "sample_preparation_id", "calibrant_run_id",
                           "calibrant_preparation_id", "source"), PairRecord),
    "scale": _object({"direction": _text, "S": _QUANTITY, "gains": _list(_number, 2), "storage_scales": _list(_text, 2)}),
}
_SIGNAL = {"name": _text, "bounds_ppm": _list(_number, 2), "evidence": _EVIDENCE, "t1_s": _number, "t1_source": _text}
_SPECTRAL_RUN = {"file": _text, "noise_bounds_ppm": _list(_number, 2), "recovery": _RECOVERY,
                 "reference": _object({"shift_ppm": _number, "source": _text})}
_CONDITIONS = {**_texts("area_mode", "probe", "pulse_sequence", "decoupling", "solvent", "tube", "instrument",
                        "excitation", "fill", "reader_scale_source"),
               **dict.fromkeys(("scans", "receiver_gain", "temperature", "reader_scale", *_RECORDED), _number)}


def _gap(name, reason="is missing"):
    return Prerequisite(name, "cannot be assessed", reason)


def _components(evidence, path, gaps):
    values = dict(evidence)
    for group in ("target_lines", "foreign_lines"):
        lines = []
        for index, line in enumerate(values.get(group, ())):
            item = f"{path}.{group}[{index}]"
            missing = [name for name in _LINE_FIELDS if line.get(name) is None]
            if missing:
                gaps.extend(_gap(f"{item}.{name}") for name in missing)
            else:
                lines.append(_build(Line, line, item))
        values[group] = tuple(lines)
    if values.setdefault("satellite_fraction", None) is None and values.get("satellite_convention") in ("include", "exclude"):
        gaps.append(_gap(f"{path}.satellite_fraction"))
    return _build(ComponentEvidence, values, path)


def _pair_id_gaps(record, pair, roles):
    ids = [record.get(role, {}).get("run_id") for role in roles]
    for role, run_id in zip(roles, ids):
        if _missing(run_id):
            yield _gap(f"{role} run_id")
        elif pair is not None:
            name = f"{role}_run_id"
            declared = getattr(pair, name)
            if declared != run_id:
                yield _gap(f"pair.{name}", "is missing" if _missing(declared) else f"must match {role}.run_id")
    if not _missing(ids[0]) and ids[0] == ids[1]:
        yield _gap("run_id", "sample and calibrant run_id must differ")


def parse_record(data: dict, base_dir: str | PathLike) -> AssayRecord:
    """Check a decoded record and build the pipeline inputs; relative paths resolve against base_dir."""
    if not isinstance(data, dict):
        _fail("", "must be an object")
    if data.get("schema") != SCHEMA:
        _fail("schema", f"must be {SCHEMA!r}")
    mode = data.get("mode")
    if mode not in MODES:
        _fail("mode", f"must be one of {', '.join(MODES)}")
    # Entered areas are calculated with the internal-standard model only.
    entered = mode == "calculation_only"
    method = "internal_standard" if entered else mode
    internal, roles = method == "internal_standard", _ROLES[method]
    run = {"run_id": _text, "preparation_id": _text, **({} if entered else _SPECTRAL_RUN)}
    if not internal:
        run["conditions"] = _object(_CONDITIONS)
    signal = {"name": _text, "area": _number, "area_source": _text} if entered else _SIGNAL
    record = _object({**_SPEC, **({} if internal else _EXTERNAL),
                      **dict.fromkeys(("run",) if internal else roles, _object(run)),
                      "quantities": _object(dict.fromkeys(_QUANTITIES[method], _QUANTITY)),
                      "signals": _object(dict.fromkeys(roles, _list(_object(signal))))})(data, "")

    signals, seen, inputs, gaps = record.get("signals", {}), {}, {}, []
    for role in roles:
        declared = []
        for index, entry in enumerate(signals.get(role, ())):
            item, name = f"signals.{role}[{index}]", entry.get("name")
            path = f"{item}.name"
            if not _has_text(name) or name in seen:
                _fail(path, f"repeats {name!r} from {seen[name]}" if name in seen else "must be nonblank text")
            seen[name] = path
            if not entered:
                evidence = entry.get("evidence")
                evidence = _components(evidence, f"{item}.evidence", gaps) if evidence else None
                declared.append(_build(DeclaredSignal, {**entry, "evidence": evidence}, item))
        if not entered:
            inputs[role] = tuple(declared)

    pair, scale = record.get("pair"), record.get("scale", {})
    if not internal:
        gaps.extend(_pair_id_gaps(record, pair, roles))
        direction = scale.get("direction")
        if direction is None:
            gaps.append(_gap("scale.direction"))
        elif direction != "sample/calibrant":
            _fail("scale.direction", "must be 'sample/calibrant'")
    base_dir = Path(base_dir).resolve()
    runs = {}
    for role in ("run",) if internal else roles:
        run = runs[role] = record.get(role, {})
        if entered:
            continue
        if not _has_text(run.get("file")):
            _fail(f"{role}.file", "must name the data file")
        run["file"] = Path(base_dir, run["file"])
        shift, source = (run.get("reference", {}).get(key) for key in ("shift_ppm", "source"))
        prefix = "" if internal else f"{role}_"
        if shift is None or _missing(source):
            gaps.append(_gap("reference" if internal else f"{role} reference", "Declare shift_ppm with its source"))
        else:
            inputs[f"{prefix}reference"] = (shift, source)
        inputs.update({f"{prefix}fid": run["file"], f"{prefix}noise_bounds_ppm": run.get("noise_bounds_ppm"),
                       f"{prefix}recovery": run.get("recovery") or RecoveryEvidence(extra_recovery_s=None),
                       **{f"{prefix}{key}": run.get(key) for key in ("run_id", "preparation_id")}})

    common = {**{key: value for key, value in record.items() if key in _PASSED}, **record.get("quantities", {}),
              "correlations": (*record.get("correlations", ()), *record.get("covariances", ())),
              "prerequisites": tuple(gaps),
              **{f"{role}_signals": tuple(entry["name"] for entry in signals.get(role, ())) for role in roles}}
    if internal:
        tube = runs["run"].get("preparation_id")
        areas = {f"{role}_{name}": tuple(entry.get(key) for entry in signals.get(role, ())) for role in roles
                 for name, key in (("areas", "area"), ("area_sources", "area_source"))} if entered else {}
        assay = _build(InternalStandardAssay, {**common, **areas, "analyte_tube_id": tube, "standard_tube_id": tube}, "")
    else:
        # The pipeline identifies a spectral run by its file, so the file stands in the assay's run_id.
        external = {role: _build(ExternalRun, {**run.get("conditions", {}), "preparation_id": run.get("preparation_id"),
                                               "run_id": str(run["file"])}, role) for role, run in runs.items()}
        if pair is not None:
            ids = {f"{role}_run_id": external[role].run_id for role, run in runs.items()
                   if not _missing(run.get("run_id")) and getattr(pair, f"{role}_run_id") == run["run_id"]}
            source = join(base_dir, pair.source) if _has_text(pair.source) else None
            pair = replace(pair, **ids, source=source)
        assay = _build(ExternalCalibrationAssay, {**common, **external, "pair": pair, "S": scale.get("S"), "S_gains":
                                                  scale.get("gains"), "S_storage_scales": scale.get("storage_scales")}, "")
    return AssayRecord(mode, assay, runs, inputs)


def load_record(path: str | PathLike) -> AssayRecord:
    """Read a record file; relative paths in it resolve against the file's folder."""
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"), parse_float=_json_float, parse_int=_json_int,
                          object_pairs_hook=_json_object)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RecordError(f"record: invalid JSON: {error}") from None
    return parse_record(data, path.parent)


def assess_record(record: AssayRecord):
    """Return the pipeline assessment for a spectral mode, or the calculation-only result, which never passes."""
    if record.mode == "calculation_only":
        return calculate_internal_standard(record.assay)
    route = assess_internal_standard if record.mode == "internal_standard" else assess_external_calibration
    return route(assay=record.assay, **record.inputs)
