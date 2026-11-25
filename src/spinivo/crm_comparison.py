"""Compare one checked assay with a laboratory's certified reference material."""

import json
from math import isclose, isfinite
from os import PathLike
from pathlib import Path

from . import __version__, accuracy, assay_record, assessment_report
from .assay_record import RecordError, _fail, _json_float, _json_int, _json_object, _number, _object, _text, _texts
from .external_calibration import _missing

__all__ = ["compare", "render_json", "render_text"]

SCHEMA = "spinivo.crm-comparison/1"
REPORT_SCHEMA = "spinivo.crm-comparison-report/1"
LIMITATION = "Neither verdict qualifies a method or releases a batch."
render_json = assessment_report.render_json


def _positive(value, path):
    value = _number(value, path)
    if value is not None and value <= 0:
        _fail(path, "must be positive")
    return value


def _mass_fraction(value, path):
    value = _positive(value, path)
    if value is not None and value > 100:
        _fail(path, "must be at most 100")
    return value


def _uncertainty(value, path):
    value = _number(value, path)
    if value is not None and value < 0:
        _fail(path, "must be nonnegative")
    return value


def _boolean(value, path):
    if value is not None and not isinstance(value, bool):
        _fail(path, "must be a boolean or null")
    return value


_THRESHOLD = _object({"value": _positive, "basis": _text})
_COMPARISON = _object({
    **_texts("schema", "comparison_id", "assay_record"),
    "crm": _object({**_texts("id", "lot", "certificate_source"), "mass_fraction_pct": _mass_fraction,
                    "standard_u_pp": _uncertainty, "expanded_u_pp": _uncertainty, "coverage_factor": _positive}),
    "independence": _object({"independent_of_calibration": _boolean, "basis": _text}),
    "comparison_k": _THRESHOLD,
    "lab_limit_pp": _THRESHOLD,
}, allow_null=False)


def _load_comparison(path):
    try:
        data = json.loads(path.read_text(encoding="utf-8"), parse_float=_json_float, parse_int=_json_int,
                          object_pairs_hook=_json_object)
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise RecordError(f"record: invalid JSON: {error}") from None
    comparison = _COMPARISON(data, "")
    if comparison.get("schema") != SCHEMA:
        _fail("schema", f"must be {SCHEMA!r}")
    return comparison


def _reference_u(crm):
    standard, expanded, factor = (crm.get(key) for key in ("standard_u_pp", "expanded_u_pp", "coverage_factor"))
    if standard is not None:
        if expanded is not None:
            _fail("crm", "supply either standard_u_pp or expanded_u_pp, not both")
        return standard, "CRM certificate standard uncertainty", None
    if expanded is None:
        return None, None, "Missing crm.standard_u_pp or crm.expanded_u_pp"
    if factor is None:
        return None, None, "Missing crm.coverage_factor for expanded_u_pp"
    standard = expanded / factor
    if not isfinite(standard) or (expanded and standard == 0):
        _fail("crm", "expanded_u_pp / coverage_factor is outside the float range")
    return standard, "CRM certificate expanded uncertainty divided by its coverage factor", None


def _checked_assay(folder, file):
    if file is None or not file.strip():
        return {}
    if Path(file).is_absolute() or "\0" in file:
        _fail("assay_record", "must be a valid path relative to the comparison file")
    try:
        record = assay_record.load_record(folder / file)
    except UnicodeEncodeError:
        _fail("assay_record", "path cannot be encoded")
    result = assay_record.assess_record(record)
    return assessment_report.assessment_report(result, record)


def _verdict(error, threshold, reasons, criterion):
    if reasons:
        return {"status": "cannot be assessed", "reason": ". ".join(reasons)}
    if error <= threshold or isclose(error, threshold, rel_tol=1e-12, abs_tol=1e-12):
        return {"status": "pass", "reason": f"Absolute error is at or below {criterion}"}
    return {"status": "fail", "reason": f"Absolute error exceeds {criterion}"}


def compare(path: str | PathLike) -> dict:
    """Run the referenced assay and return the versioned CRM comparison report."""
    path = Path(path)
    data = _load_comparison(path)
    crm, independence = data.get("crm", {}), data.get("independence", {})
    reference_u, reference_basis, reference_gap = _reference_u(crm)
    assay = _checked_assay(path.parent, data.get("assay_record"))
    status = assay.get("status")
    purity = assay["purity"]["value"] if status == "pass" else None
    standard_u = assay["method_report"]["u_c"] if status == "pass" else None
    measured = purity * 100 if purity is not None and isfinite(purity * 100) else None
    measurement_u = standard_u * 100 if standard_u is not None and isfinite(standard_u * 100) else None

    comparison_id, reference = data.get("comparison_id"), crm.get("mass_fraction_pct")
    independent, independence_basis = independence.get("independent_of_calibration"), independence.get("basis")
    k, k_basis = (data.get("comparison_k", {}).get(key) for key in ("value", "basis"))
    limit, limit_basis = (data.get("lab_limit_pp", {}).get(key) for key in ("value", "basis"))
    common = []
    if not assay:
        common.append("Missing assay_record")
    elif status != "pass":
        common.append(f"Assay status is {status}. A checked pass is required")
    elif measured is None:
        common.append("Checked purity is unavailable in percent")
    if reference is None:
        common.append("Missing crm.mass_fraction_pct")
    if independent is not True:
        common.append("independence.independent_of_calibration must be declared true")
    if _missing(independence_basis):
        common.append("Missing independence.basis")

    result = {}
    if measured is not None and reference is not None:
        # A signed difference remains descriptive even without evidence of independence.
        identifier = comparison_id if not _missing(comparison_id) else "comparison"
        preparation = accuracy.Preparation(identifier, assay["method"], "unspecified", "unspecified")
        row = accuracy.AccuracyRow(identifier, preparation, "pass", identifier, measured, reference,
                                   measurement_u if not common else None, reference_u, reference_basis)
        # Use an explicit placeholder tolerance if no lab limit was supplied. Only the row metrics are used.
        scored = accuracy.summarize_accuracy([row], preparations=[preparation], tolerance_pp=limit or 0.0, k=k)
        result = scored["preparations"][0]

    error, error_u = result.get("error_pp"), result.get("error_standard_u_pp")
    absolute = abs(error) if error is not None else None
    uncertainty_underflow = error_u == 0 and (measurement_u != 0 or reference_u != 0)
    if uncertainty_underflow:
        error_u = None
    threshold = k * error_u if k is not None and error_u is not None else None
    agreement_gaps = [reason for missing, reason in (
        (k is None, "Missing comparison_k.value"), (_missing(k_basis), "Missing comparison_k.basis"),
        (reference_gap is not None, reference_gap), (measurement_u is None, "Missing measurement standard uncertainty"),
    ) if missing]
    if uncertainty_underflow:
        agreement_gaps.append("Signed-error variance is below the float range")
    elif not common and not agreement_gaps and error_u is None:
        agreement_gaps.append(result["error_u_reason"])
    if threshold is not None and not isfinite(threshold):
        threshold = None
        agreement_gaps.append("k u(error) is outside the float range")
    limit_gaps = [reason for missing, reason in (
        (limit is None, "Missing lab_limit_pp.value"), (_missing(limit_basis), "Missing lab_limit_pp.basis"),
    ) if missing]

    ratio = result.get("standard_u_ratio")
    return {
        "schema": REPORT_SCHEMA, "spinivo_version": __version__, "comparison_id": comparison_id,
        "assay": {"assay_id": assay.get("assay_id"), "method": assay.get("method"),
                  "status": status, "reasons": assay.get("reasons", [])},
        "crm": {key: crm.get(key) for key in ("id", "lot", "certificate_source")},
        "independence": {"independent_of_calibration": independent, "basis": independence_basis},
        "measured_purity_pct": measured, "reference_purity_pct": reference,
        "measurement_standard_u_pp": measurement_u, "reference_standard_u_pp": reference_u,
        "error_pp": error, "absolute_error_pp": absolute, "error_standard_u_pp": error_u,
        "k": k, "k_basis": k_basis, "comparison_u_pp": threshold,
        "standard_u_ratio": ratio if threshold and ratio is not None and isfinite(ratio) else None,
        "lab_limit_pp": limit, "lab_limit_basis": limit_basis,
        "uncertainty_agreement": _verdict(absolute, threshold, common or agreement_gaps, "k u(error)"),
        "lab_limit": _verdict(absolute, limit, common or limit_gaps, "the laboratory limit"),
        "limitation": LIMITATION,
    }


def _escape_controls(message):
    return "".join(char if char.isprintable() else ascii(char)[1:-1] for char in message)


def render_text(report: dict) -> str:
    """Render the same comparison values and verdicts in plain text."""
    assay, crm, independence = (report[key] for key in ("assay", "crm", "independence"))
    lines = []
    for label, value, unit in (
        ("Comparison", report["comparison_id"], ""),
        ("Report", f'{report["schema"]}, spinivo {report["spinivo_version"]}', ""),
        ("Assay", assay["assay_id"], ""), ("Assay method", assay["method"], ""), ("Assay status", assay["status"], ""),
        ("CRM", crm["id"], ""), ("CRM lot", crm["lot"], ""), ("Certificate source", crm["certificate_source"], ""),
        ("Independent of calibration", independence["independent_of_calibration"], ""),
        ("Independence basis", independence["basis"], ""),
        ("Measured mass fraction", report["measured_purity_pct"], " %"),
        ("Reference mass fraction", report["reference_purity_pct"], " %"),
        ("Measurement standard uncertainty", report["measurement_standard_u_pp"], " pp"),
        ("Reference standard uncertainty", report["reference_standard_u_pp"], " pp"),
        ("Signed error", report["error_pp"], " pp"), ("Absolute error", report["absolute_error_pp"], " pp"),
        ("u(error)", report["error_standard_u_pp"], " pp"), ("k", report["k"], ""), ("k basis", report["k_basis"], ""),
        ("k u(error)", report["comparison_u_pp"], " pp"), ("|error| / (k u(error))", report["standard_u_ratio"], ""),
        ("Lab limit", report["lab_limit_pp"], " pp"), ("Lab limit basis", report["lab_limit_basis"], ""),
    ):
        lines.append(f"{label}: " + ("missing" if value is None else f"{value}{unit}"))
    lines.extend(f"Assay reason: {reason}" for reason in assay["reasons"])
    for key, label in (("uncertainty_agreement", "Uncertainty agreement"), ("lab_limit", "Lab limit verdict")):
        verdict = report[key]
        lines.append(f'{label}: {verdict["status"]}. {verdict["reason"]}')
    lines.append(report["limitation"])
    return "\n".join(_escape_controls(line) for line in lines)
