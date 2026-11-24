"""Render a same-tube purity result as versioned JSON and plain text."""

from decimal import Decimal

from .integral_report import _json_value, render_json
from .internal_standard import Correlation, Covariance, InternalStandardAssay, InternalStandardResult, _finite, _QUANTITIES
from .pipeline import InternalStandardAssessment

__all__ = ["purity_report", "render_json", "render_text"]


def _text_or_none(value):
    return value if isinstance(value, str) else None


def _input(value, u, unit, source, source_ids=()):
    record = {"value": value if _finite(value) else None, "u": u if _finite(u) else None,
              "unit": _text_or_none(unit), "source": _text_or_none(source)}
    if isinstance(source_ids, (tuple, list)) and source_ids:
        record["source_ids"] = [_text_or_none(item) for item in source_ids]
    return record


def _covariances(correlations, budget):
    rows = []
    for index, c in enumerate(correlations if isinstance(correlations, (tuple, list)) else ()):
        if not isinstance(c, (Correlation, Covariance)):
            continue
        pair = [_text_or_none(name) for name in c.inputs] if isinstance(c.inputs, (tuple, list)) else None
        row = budget.covariances[index] if budget else None
        r = c.r if isinstance(c, Correlation) else row.r if row else None
        evidence = {"inputs": pair, "r": r if _finite(r) else None, "source": _text_or_none(c.source),
                    "term": row.term if row else None, "share": row.share if row else None}
        if isinstance(c, Covariance):
            evidence["covariance"] = c.value if _finite(c.value) else None
            evidence["units"] = [_text_or_none(unit) for unit in c.units] if isinstance(c.units, (tuple, list)) else None
        if c.independence_reason is not None:
            evidence["independence_reason"] = _text_or_none(c.independence_reason)
        rows.append(evidence)
    return rows


def purity_report(result: InternalStandardAssessment | InternalStandardResult,
                  assay: InternalStandardAssay | None = None) -> dict:
    """Render a pipeline assessment or a calculation-only result with its record."""
    if isinstance(result, InternalStandardAssessment):
        if assay is not None:
            raise TypeError("Pass the assessment by itself")
        assay, result = result.assay, result.result
    elif result.status == "pass":
        raise ValueError("A checked pass report requires an InternalStandardAssessment")
    elif assay is None:
        raise TypeError("A calculation result requires its assay")
    budget = result.budget
    relative_u = assay.R_relative_u
    ratio_u = abs(result.R) * relative_u if _finite(result.R) and _finite(relative_u) and relative_u >= 0 else None
    inputs = {"R": _input(result.R, ratio_u, "ratio", assay.R_u_source, assay.R_source_ids),
              "R_relative_u": _input(relative_u, None, "fraction", assay.R_u_source),
              "N_a": _input(assay.N_a, None, "count", None),
              "N_s": _input(assay.N_s, None, "count", None)}
    for role in ("analyte", "standard"):
        name = f"{role}_areas"
        areas = getattr(assay, name)
        sources = getattr(assay, f"{role}_area_sources")
        aligned_sources = None
        if isinstance(sources, (tuple, list)):
            area_count = len(areas or ())
            aligned_sources = [_text_or_none(source) for source in sources[:area_count]]
            aligned_sources.extend([None] * (area_count - len(aligned_sources)))
        inputs[name] = {
            "value": [area if _finite(area) else None for area in areas] if areas is not None else None,
            "u": None, "unit": _text_or_none(assay.area_unit), "source": aligned_sources,
        }
    for name in _QUANTITIES:
        q = getattr(assay, name)
        inputs[name] = _input(q.value, q.standard_u, q.unit, q.source, q.source_ids) if q else _input(None, None, None, None)

    rows = [{"name": row.name, "sensitivity": row.sensitivity, "contribution": row.contribution, "share": row.share}
            for row in budget.inputs] if budget else []
    return _json_value({
        "schema": "spinivo.internal_standard_purity/2",
        "method": "internal standard, same tube",
        "assay_id": result.assay_id, "analyte_tube_id": assay.analyte_tube_id, "standard_tube_id": assay.standard_tube_id,
        "status": result.status, "purity": {"value": result.value, "unit": "fraction"},
        "I_a": result.I_a, "I_s": result.I_s, "R": result.R, "inputs": inputs,
        "budget": rows, "covariances": _covariances(assay.correlations, budget), "variance": budget.variance if budget else None,
        "u_c": budget.standard_u if budget else None, "u_c/w": budget.relative_standard_u if budget else None,
        "k": {"value": assay.k.value if _finite(assay.k.value) else None,
              "source": _text_or_none(assay.k.source)} if assay.k else None,
        "U": budget.expanded_u if budget else None,
        "U_reason": budget.expanded_u_reason if budget else "purity cannot be assessed",
        "prerequisites": assay.prerequisites, "reasons": result.reasons,
    })


def _decimal_places(value):
    if not value:
        return 2
    rounded = round(value, 1 - value.adjusted())
    return 1 - rounded.adjusted()


def _percent_text(percent, places):
    if places >= 0:
        return f"{percent:.{places}f} %"
    rounded = Decimal(f"{percent.scaleb(places):.0f}").scaleb(-places)
    return f"{rounded:f} %"


def _number_text(value):
    return f"{value:.8g}" if value is not None else "unavailable"


def render_text(report: dict) -> str:
    """Show calculated purity in percent after the status reasons and keep independent areas visible."""
    lines = [f'{report["assay_id"]} ({report["method"]}): {report["status"]}']
    lines.extend(f"  {name} = {_number_text(report[name])}" for name in ("I_a", "I_s", "R"))
    lines.extend(f"  {reason}" for reason in report["reasons"])
    if report["purity"]["value"] is None:
        return "\n".join(lines)

    standard = Decimal(str(report["u_c"])).scaleb(2)
    purity = Decimal(str(report["purity"]["value"])).scaleb(2)
    places = _decimal_places(standard)
    label = "" if report["status"] == "pass" else " (calculation only, not checked)"
    lines.append(f"  Purity = {_percent_text(purity, places)}, u_c = {_percent_text(standard, places)}{label}")
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
