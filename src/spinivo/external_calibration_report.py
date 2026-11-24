"""Render an external-calibration purity result as versioned JSON and plain text."""

from dataclasses import fields
from decimal import Decimal

from .external_calibration import (
    ExternalCalibrationAssay, ExternalCalibrationResult, PairRecord, _PREPARATION, _RECORDED,
    _missing, _scan_area,
)
from .integral_report import _json_value, render_json
from .internal_standard import Quantity, _finite
from .internal_standard_report import (
    _covariances, _decimal_places, _input, _number_text, _percent_text, _text_or_none,
)
from .pipeline import ExternalCalibrationAssessment

__all__ = ["purity_report", "render_json", "render_text"]


def _run_report(assay, role, ids=None):
    """Keep scan-normalized measured areas even when reader correction is unavailable.

    An assessment's run IDs replace the entered run_id, which the record route sets to the file path.
    """
    run = getattr(assay, role)
    numeric = ("area", "scans", "receiver_gain", "temperature", "reader_scale", "calibrated_flip_angle_deg", *_RECORDED)
    record = {}
    for field in fields(run):
        value = getattr(run, field.name)
        if field.name in numeric:
            record[field.name] = value if _finite(value) else None
        else:
            record[field.name] = None if _missing(value) else _text_or_none(value)
    record["per_scan_area"] = _scan_area(run)
    header = {name: value for owner, name, value in assay._header_values if owner == role}
    record["header"] = header or None
    record["entered"] = None
    for owner, entered, _ in assay._linked_runs:
        if owner == role:
            values = {**vars(entered), **({} if ids is None else {"run_id": ids.run_id})}
            record["entered"] = {name: value for name, value in values.items()
                                 if _finite(value) or isinstance(value, str) and not _missing(value)} or None
    return record


def _pair_report(pair):
    if not isinstance(pair, PairRecord):
        return None
    return {name: None if _missing(value) else _text_or_none(value) for name, value in vars(pair).items()}


def purity_report(result: ExternalCalibrationAssessment | ExternalCalibrationResult,
                  assay: ExternalCalibrationAssay | None = None) -> dict:
    """Report the calculation's record. The optional assay is accepted for existing callers."""
    runs = {}
    if isinstance(result, ExternalCalibrationAssessment):
        runs = {run.run: run for run in result.runs}
        result = result.result  # Its calculation record is the assessment's linked record.
    elif result.status == "pass":
        raise ValueError("A checked pass report requires an ExternalCalibrationAssessment")
    if result._assay is None:
        raise ValueError("External calibration report requires the calculation record")
    assay = result._assay
    budget = result.budget
    relative_u = assay.R_relative_u
    ratio_u = abs(result.R) * relative_u if _finite(result.R) and _finite(relative_u) and relative_u >= 0 else None
    inputs = {"R": _input(result.R, ratio_u, "ratio", assay.R_u_source, assay.R_source_ids),
              "R_relative_u": _input(relative_u, None, "fraction", assay.R_u_source),
              "N_a": _input(assay.N_a, None, "count", None),
              "N_s": _input(assay.N_s, None, "count", None)}
    path = ("c_s",) if assay.c_s_route == "certified" or assay.c_s is not None else ()
    if assay.c_s_route != "certified" or any(getattr(assay, name) is not None for name in _PREPARATION):
        path += _PREPARATION
    for name in (*path, "M_a", "V_a", "m_a", "K", "D", "S"):
        q = getattr(assay, name)
        inputs[name] = (_input(q.value, q.standard_u, q.unit, q.source, q.source_ids)
                        if isinstance(q, Quantity) else _input(None, None, None, None))

    rows = [{"name": row.name, "sensitivity": row.sensitivity, "contribution": row.contribution, "share": row.share}
            for row in budget.inputs] if budget else []
    return _json_value({
        "schema": "spinivo.external_calibration_purity/3",
        "method": "external calibration",
        "assay_id": _text_or_none(result.assay_id),
        "status": result.status, "purity": {"value": result.value, "unit": "fraction"},
        "J_a": result.J_a, "J_s": result.J_s, "R": result.R,
        "c_s_route": _text_or_none(assay.c_s_route),
        "c_s": result.c_s, "c_a": result.c_a, "concentration_unit": result.concentration_unit,
        **{role: _run_report(assay, role, runs.get(role)) for role in ("sample", "calibrant")},
        "pair": _pair_report(assay.pair),
        "S_gains": ([gain if _finite(gain) else None for gain in assay.S_gains]
                    if isinstance(assay.S_gains, (tuple, list)) else None),
        "S_storage_scales": ([_text_or_none(scale) for scale in assay.S_storage_scales]
                             if isinstance(assay.S_storage_scales, (tuple, list)) else None),
        "tolerances": dict(result.transfer_tolerances),
        "K_covers": ([_text_or_none(name) for name in assay.K_covers]
                     if isinstance(assay.K_covers, (tuple, list)) else None),
        "inputs": inputs, "budget": rows, "covariances": _covariances(assay.correlations, budget),
        "variance": budget.variance if budget else None,
        "u_c": budget.standard_u if budget else None, "u_c/w": budget.relative_standard_u if budget else None,
        "k": {"value": assay.k.value if _finite(assay.k.value) else None,
              "source": _text_or_none(assay.k.source)} if assay.k else None,
        "U": budget.expanded_u if budget else None,
        "U_reason": budget.expanded_u_reason if budget else "purity cannot be assessed",
        "prerequisites": assay.prerequisites, "reasons": result.reasons,
    })


def render_text(report: dict) -> str:
    """Show calculated purity in percent after the status reasons and keep per-scan areas visible."""
    lines = [f'{report["assay_id"]} ({report["method"]}): {report["status"]}']
    for role in ("sample", "calibrant"):
        run = report[role]
        area, scans = (_number_text(run[name]) for name in ("area", "scans"))
        per_scan = f', area per scan = {_number_text(run["per_scan_area"])}' if "per_scan_area" in run else ""
        lines.append(f"  {role} area = {area}, scans = {scans}{per_scan}")
    lines.extend(f"  {name} = {_number_text(report[name])}" for name in ("J_a", "J_s", "R"))
    unit = report["concentration_unit"]
    if report["c_s"] is not None:
        route = f' ({report["c_s_route"]})' if report.get("c_s_route") else ""
        lines.append(f'  c_s = {_number_text(report["c_s"])} {unit}{route}')
    lines.extend(f"  {reason}" for reason in report["reasons"])
    if report["purity"]["value"] is None:
        return "\n".join(lines)

    lines.append(f'  c_a = {_number_text(report["c_a"])} {unit}')
    for name in (name for name in ("K", "S") if name in report["inputs"]):
        factor = report["inputs"][name]
        lines.append(f'  {name} = {_number_text(factor["value"])} ({factor["source"]})')

    standard = Decimal(str(report["u_c"])).scaleb(2)
    purity = Decimal(str(report["purity"]["value"])).scaleb(2)
    places = _decimal_places(standard)
    lines.append(f"  Purity = {_percent_text(purity, places)}, u_c = {_percent_text(standard, places)}")
    if report["U"] is None:
        lines.append(f'  U unavailable ({report["U_reason"]})')
    else:
        expanded = Decimal(str(report["U"])).scaleb(2)
        lines.append(f'  U = {_percent_text(expanded, _decimal_places(expanded))} with k = {report["k"]["value"]:g}')
    lines.append(f'  {"Input":<6} {"Value":>12} {"u":>12} {"Sensitivity":>12} {"Contribution":>12} {"Share":>12}')
    for row in report["budget"]:
        name = row["name"]
        quantity = report["inputs"][name]
        columns = (quantity["value"], quantity["u"], row["sensitivity"], row["contribution"], row["share"])
        lines.append(f"  {name:<6} " + " ".join(f"{_number_text(value):>12}" for value in columns))
    if report["covariances"]:
        lines.append(f'  {"Covariance":<12} {"r":>12} {"Variance term":>14} {"Share":>12}')
        for row in report["covariances"]:
            pair = ", ".join(row["inputs"])
            lines.append(f'  {pair:<12} {_number_text(row["r"]):>12} '
                         f'{_number_text(row["term"]):>14} {_number_text(row["share"]):>12}')
    return "\n".join(lines)
