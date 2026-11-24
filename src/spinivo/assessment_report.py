"""Render what assess_record returns as one versioned JSON report (spinivo.assessment/1) and plain text."""

from . import __version__, external_calibration_report, internal_standard_report
from .assay_record import AssayRecord
from .integral_checks import IntegralThresholds
from .integral_report import _json_value, _observation_report, region_report, render_json, render_text as region_text
from .internal_standard import InternalStandardResult, calculate_internal_standard
from .internal_standard_report import _text_or_none
from .pipeline import ExternalCalibrationAssessment, InternalStandardAssessment

__all__ = ["assessment_report", "render_json", "render_text"]

SCHEMA = "spinivo.assessment/1"
NOT_ASSESSED = ("raw data, processing and referencing", "acquisition checks: recovery, snr, resolution",
                "integral region checks", "capture pair checks")
_METHODS = {"internal_standard": internal_standard_report, "external_calibration": external_calibration_report}
_UNITS = (("_mhz", "MHz"), ("_hz", "Hz"), ("_ppm", "ppm"), ("_us", "us"), ("_deg", "deg"), ("_s", "s"),
          ("_points", "count"), ("scans", "count"), ("height", "intensity"), ("noise_sd", "intensity"))


def _check_report(check):
    inputs = [{"name": item.name, "value": item.value, "source": item.source,
               "unit": next((unit for suffix, unit in _UNITS if item.name.endswith(suffix)), None)}
              for item in check.inputs]
    return {"run": check.run, "signal": check.signal, "check": check.check, "status": check.status,
            "value": check.value, "unit": check.unit, "limit": check.limit, "reason": check.reason, "inputs": inputs}


def _region_entry(source, role, signal, region, profiles):
    if region is not None and region.baseline.details.recovery.limit != profiles.baseline_limit:
        raise ValueError("Baseline limit does not match the region observations")
    report = None if region is None else region_report(region, run=source, signal=signal, thresholds=profiles.integral)
    return {"role": role, "signal": signal, "report": report}


def _run_report(run, profiles):
    acquisition, regions = run.acquisition, run.regions
    return {"source": run.source, "run_id": _text_or_none(run.run_id),
            "preparation_id": _text_or_none(run.preparation_id), "processing": run.settings,
            "acquisition": None if acquisition is None else [_check_report(check) for check in acquisition],
            "regions": None if regions is None else [_region_entry(run.source, *entry, profiles) for entry in regions]}


def assessment_report(obj, record: AssayRecord | None = None, *, thresholds: IntegralThresholds | None = None) -> dict:
    """Report with the profiles stored at evaluation; a calculation-only result needs its record.

    A threshold override that differs from the applied profile raises. Purity is a result only on pass.
    """
    if isinstance(obj, (InternalStandardAssessment, ExternalCalibrationAssessment)):
        method = mode = "internal_standard" if isinstance(obj, InternalStandardAssessment) else "external_calibration"
        method_report, profiles, runs = _METHODS[method].purity_report(obj), obj.profiles, obj.runs
        pairs = [{"signals": (first, second), **_observation_report(pair)} for first, second, pair in obj.capture_pairs]
    elif isinstance(obj, InternalStandardResult):
        if record is None:
            raise TypeError("A calculation-only result requires its record")
        if record.mode != "calculation_only":
            raise ValueError(f"Record mode {record.mode!r} does not match the calculation_only assessment")
        rebuilt = calculate_internal_standard(record.assay)
        if rebuilt != obj:
            raise ValueError("Calculation result does not match the supplied record")
        method, mode, profiles, runs = "internal_standard", "calculation_only", None, ()
        method_report, pairs = internal_standard_report.purity_report(rebuilt, record.assay), None
    else:
        raise TypeError(f"Unsupported assessment: {type(obj).__name__}")
    if record is not None and record.mode != mode:
        raise ValueError(f"Record mode {record.mode!r} does not match the {mode} assessment")
    if thresholds is not None and (profiles is None or thresholds not in (profiles.integral, profiles.requested_integral)):
        raise ValueError("Threshold override does not match the profile applied at evaluation")
    status = method_report["status"]
    return _json_value({
        "schema": SCHEMA, "spinivo_version": __version__, "method": method, "mode": mode,
        "assay_id": method_report["assay_id"], "status": status,
        "purity": {"value": method_report["purity"]["value"] if status == "pass" else None, "unit": "fraction"},
        "reasons": method_report["reasons"],
        "not_assessed": NOT_ASSESSED if profiles is None else (),
        "profiles": None if profiles is None else {
            "integral": profiles.integral, "acquisition": profiles.acquisition,
            "baseline": {"algorithm_version": profiles.baseline_algorithm, "limit": profiles.baseline_limit},
        },
        "runs": None if profiles is None else {run.run: _run_report(run, profiles) for run in runs},
        "capture_pairs": pairs,
        "method_report": method_report,
        **({"baseline_shape": obj.baseline_shape} if profiles is not None else {}),
    })


def _text(value):
    if isinstance(value, list):
        return f'[{", ".join(map(_text, value))}]'
    return "missing" if value is None else f"{value:.8g}" if isinstance(value, float) else str(value)


def _pairs(values):
    return ", ".join(f"{key} = {_text(value)}" for key, value in values.items())


def _check_text(name, check):
    return (f'{name}: {check["status"]}, value = {_text(check["value"])} {check["unit"]}, '
            f'limit = {_text(check["limit"])}: {check["reason"]}')


def baseline_shape_text(shape):
    value = shape["relative_standard_uncertainty_percent"]
    detail = f"{_text(value)} % relative standard uncertainty" if value is not None else shape["reason"]
    return f'Baseline shape ({shape["method_version"]}, diagnostic): {shape["status"]}, {detail}'


def _run_text(label, run):
    lines = [f'  Run {label}: {_text(run["source"])}, run_id = {_text(run.get("run_id"))}, '
             f'preparation_id = {_text(run.get("preparation_id"))}',
             f'    processing: {_pairs(run["processing"]) if run["processing"] else "missing"}',
             *(f"    {part}: not measured" for part in ("acquisition", "regions") if run[part] is None)]
    for check in run["acquisition"] or ():
        lines.append("    " + _check_text(f'{check["signal"]} {check["check"]}', check))
        lines.extend(f'      {item["name"]} = {_text(item["value"])}'
                     f'{" " + item["unit"] if item["unit"] and item["value"] is not None else ""} '
                     f'({_text(item["source"])})' for item in check["inputs"])
    for region in run["regions"] or ():
        if region["report"] is None:
            lines.append(f'    {region["signal"]} {region["role"]}: no region result')
            continue
        lines.extend(f"    {line}" for line in region_text(region["report"]).splitlines())
        lines.extend("      " + _check_text(obs["check"], obs) for obs in region["report"]["observations"])
    return lines


def render_text(report: dict) -> str:
    """Keep every pass, flag and gap record, the blocking reasons and what the mode did not assess."""
    purity, method = report["purity"], report["method_report"]
    checked = f'{_text(purity["value"])} {purity["unit"]}' if purity["value"] is not None else report["status"]
    lines = [f'{_text(report["assay_id"])} ({report["method"]}, {report["mode"]} mode): {report["status"]}',
             f'  {report["schema"]}, spinivo {report["spinivo_version"]}', f"  Checked purity: {checked}"]
    if report.get("baseline_shape") is not None:
        lines.append("  " + baseline_shape_text(report["baseline_shape"]))
    for title, items in (("Blocking reasons", report["reasons"]),
                         (f'Not assessed in {report["mode"]} mode', report["not_assessed"]),
                         ("Checks", [f'{c["name"]}: {c["status"]}: {c["reason"]}' for c in method["prerequisites"]])):
        if items:
            lines.extend((f"  {title}", *(f"    {item}" for item in items)))
    lines.extend(f"  {name} profile: {_pairs(profile)}" for name, profile in (report["profiles"] or {}).items())
    for label, run in (report["runs"] or {}).items():
        lines.extend(_run_text(label, run))
    for pair in report.get("capture_pairs") or ():
        name = "/".join(pair["signals"]) + " capture_pair"
        lines.append(f'  {_check_text(name, pair)} (threshold version = {pair["threshold_version"]})')
    lines.extend(f"  {line}" for line in _METHODS[report["method"]].render_text(method).splitlines())
    return "\n".join(lines)
