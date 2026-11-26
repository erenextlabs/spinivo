"""Render a spinivo.assessment/1 report as a self-contained HTML page."""

import os
import posixpath
import re
from decimal import Decimal
from html import escape
from pathlib import Path

from .assessment_report import _check_text, _pairs, _text, baseline_shape_text
from .integral_report import render_text as region_text
from .internal_standard_report import _decimal_places, _number_text, _percent_text

__all__ = ["render_html"]

_STYLE = """
body {
  margin: 0 auto;
  max-width: 72rem;
  padding: 1rem clamp(1rem, 4vw, 3rem) 2rem;
  font: 1rem/1.5 system-ui, sans-serif;
  color: #182630;
  background: #fff;
  overflow-wrap: anywhere;
}
h1, h2, h3, h4 { line-height: 1.2; }
h1 { margin: .25em 0; }
h2 { margin-top: 1.6em; }
header, footer { border-top: 2px solid #36686a; padding: .75rem 0; }
footer { margin-top: 2rem; font-size: .9rem; }
.purity { background: #edf5f4; padding: .75rem 1rem; font-weight: bold; }
.region-details { white-space: pre-wrap; }
.table-wrap { overflow-x: auto; margin: .75rem 0; }
table { width: 100%; min-width: 34rem; border-collapse: collapse; }
th, td { padding: .4rem .65rem; border: 1px solid #c5d0d4; text-align: left; vertical-align: top; }
thead { background: #edf2f4; }
@media print {
  body { max-width: none; padding: 0; font-size: 10pt; color: #000; }
  h1, h2, h3, h4 { break-after: avoid; }
  tr { page-break-inside: avoid; break-inside: avoid; }
  .table-wrap { overflow: visible; }
  table { min-width: 0; }
}
"""


def _paragraph(value):
    return f"<p>{escape(value)}</p>"


def _list(title, items):
    return [f"<h2>{escape(title)}</h2>", "<ul>", *(f"<li>{escape(item)}</li>" for item in items), "</ul>"]


def _table(label, headings, rows):
    header = "".join(f'<th scope="col">{escape(heading)}</th>' for heading in headings)
    return [f'<div class="table-wrap" tabindex="0" role="region" aria-label="{escape(label)}">',
            f"<table><thead><tr>{header}</tr></thead><tbody>",
            *("<tr>" + "".join(f"<td>{escape(_text(cell))}</td>" for cell in row) + "</tr>" for row in rows),
            "</tbody></table></div>"]


def _check_row(check):
    return (check["check"], check["status"], check["value"], check["unit"], check["limit"], check["reason"])


_CHECK_COLUMNS = ("Check", "Status", "Value", "Unit", "Limit", "Reason")


def _run(label, run):
    lines = [f"<h3>Run {escape(label)}</h3>",
             _paragraph(f'Source: {_text(run["source"])}. Run ID: {_text(run.get("run_id"))}. '
                        f'Preparation ID: {_text(run.get("preparation_id"))}.'),
             _paragraph(f'Processing: {_pairs(run["processing"]) if run["processing"] else "missing"}'),
             "<h4>Acquisition checks</h4>"]
    if run["acquisition"] is None:
        lines.append(_paragraph("Acquisition: not measured"))
    for check in run["acquisition"] or ():
        name = f'{check["signal"]} {check["check"]}'
        lines.append(_paragraph(_check_text(name, check)))
        if check["inputs"]:
            lines.extend(_table(f"{name} inputs", ("Input", "Value", "Unit", "Source"),
                                ((item["name"], item["value"], item["unit"], item["source"])
                                 for item in check["inputs"])))
    lines.append("<h4>Regions</h4>")
    if run["regions"] is None:
        lines.append(_paragraph("Regions: not measured"))
    for entry in run["regions"] or ():
        region = entry["report"]
        if region is None:
            lines.append(_paragraph(f'{entry["signal"]} {entry["role"]}: no region result'))
            continue
        summary, *details = region_text(region).splitlines()
        lines.append(f"<p><strong>{escape(summary)}</strong></p>")
        if details:
            details_text = "\n".join(details)
            lines.append(f'<p class="region-details">{escape(details_text)}</p>')
        lines.extend(_table(f'{entry["signal"]} observations', _CHECK_COLUMNS,
                            (_check_row(obs) for obs in region["observations"])))
    return lines


def _method_measurements(method, external):
    lines = []
    if external:
        for role in ("sample", "calibrant"):
            run = method[role]
            area, scans, per_scan = (_number_text(run[name]) for name in ("area", "scans", "per_scan_area"))
            lines.append(_paragraph(f"{role} area = {area}, scans = {scans}, area per scan = {per_scan}"))
    areas = ("J_a", "J_s", "R") if external else ("I_a", "I_s", "R")
    lines.append(_paragraph(", ".join(f"{name} = {_number_text(method[name])}" for name in areas)))
    if external:
        for name in ("c_s", "c_a"):
            if method[name] is not None:
                route = f' ({method["c_s_route"]})' if name == "c_s" and method["c_s_route"] else ""
                lines.append(_paragraph(f'{name} = {_number_text(method[name])} {method["concentration_unit"]}{route}'))
    return lines


def _purity_and_standard_text(method):
    standard = Decimal(str(method["u_c"])).scaleb(2)
    purity = Decimal(str(method["purity"]["value"])).scaleb(2)
    places = _decimal_places(standard)
    return _percent_text(purity, places), _percent_text(standard, places)


def _expanded_text(value):
    if value is None:
        return "U unavailable"
    expanded = Decimal(str(value)).scaleb(2)
    return f"U = {_percent_text(expanded, _decimal_places(expanded))}"


def _uncertainty(method, passed):
    lines = ["<h3>Uncertainty</h3>"]
    if method["purity"]["value"] is not None:
        purity, standard = _purity_and_standard_text(method)
        label = "" if passed else " (calculation only, not checked)"
        lines.append(_paragraph(f"Purity = {purity}{label}"))
        lines.append(_paragraph(f"u_c = {standard}"))
    else:
        lines.append(_paragraph("u_c unavailable"))
    expanded = _expanded_text(method["U"])
    if method["U"] is None:
        expanded += f' ({method["U_reason"]})'
    lines.append(_paragraph(expanded))
    k = method["k"]
    lines.append(_paragraph(f'k = {k["value"]:g}' if k and k["value"] is not None else "k unavailable"))
    return lines


def _input_value(value):
    return _text(value) if isinstance(value, list) else _number_text(value)


def _method_report(method, external, passed):
    lines = ["<h2>Method report</h2>", _paragraph(method["method"]), *_method_measurements(method, external),
             "<h3>Inputs</h3>"]
    lines.extend(_table("Method inputs", ("Input", "Value", "u", "Unit", "Source"),
                        ((name, _input_value(item["value"]), _input_value(item["u"]), item["unit"], item["source"])
                         for name, item in method["inputs"].items())))
    if method["budget"]:
        lines.append("<h3>Uncertainty budget</h3>")
        lines.extend(_table("Uncertainty budget", ("Input", "Sensitivity", "Contribution", "Share"),
                            ((row["name"], *(_number_text(row[key]) for key in ("sensitivity", "contribution", "share")))
                             for row in method["budget"])))
    if method["covariances"]:
        lines.append("<h3>Covariances</h3>")
        lines.extend(_table("Covariances", ("Inputs", "r", "Variance term", "Share", "Source"),
                            ((", ".join(map(_text, row["inputs"])) if row["inputs"] else None,
                              *(_number_text(row[key]) for key in ("r", "term", "share")), row["source"])
                             for row in method["covariances"])))
    lines.extend(_uncertainty(method, passed))
    return lines


def _path_spellings(path):
    quoted = repr(path)
    yield path
    yield quoted[1:-1]
    if quoted.startswith('"'):
        # A longer path containing both quotes would use single quotes in repr.
        yield quoted[1:-1].replace("'", r"\'")


def _source_names(report):
    names = {}
    for run in (report["runs"] or {}).values():
        source = run["source"]
        path = source.replace("\\", "/") if source else ""
        if "/" not in path or not path.strip("/."):
            continue
        name = posixpath.basename(posixpath.normpath(path))
        if name in ("", ".", ".."):
            name = "source"
        names.update(dict.fromkeys(_path_spellings(source), name))
    return names


def _absolute_paths(folder):
    return {os.path.abspath(folder), os.path.realpath(folder)}


def _folder_names(hide_folders):
    home = Path.home()
    if hide_folders is None:
        hide_folders = (Path.cwd(), home)
    home_paths = _absolute_paths(home)
    names = {}
    for folder in hide_folders:
        paths = _absolute_paths(folder)
        for path in paths:
            if os.path.dirname(path) != path:
                marker = "~" if path in home_paths else "."
                for spelling in {path, path.replace("\\", "/")}:
                    names.update(dict.fromkeys(_path_spellings(spelling), marker))
    return names


def _display_report(report, hide_folders):
    """Prepare display strings before formatting and escaping the HTML."""
    folders = _folder_names(hide_folders)
    names = _source_names(report) | folders
    if not names:
        return report

    # Match the longest path at each position so a shorter prefix cannot expose a nested folder.
    paths = sorted(names, key=len, reverse=True)
    path_end = r"(?![\w-]|\.[\w-])"
    source_start = r"(?<![\w./\\-])"
    folder_start = r"(?<![\w.-])"  # keeps units such as g/mol when a folder is /mol

    def start(path):
        if path not in folders:
            return source_start
        # Only a one-part folder such as /mol needs a start boundary; longer folders are hidden anywhere.
        return folder_start if not re.search(r"[/\\]", path.strip("/\\")) else ""

    pattern = re.compile("|".join(start(path) + re.escape(path) + path_end for path in paths))

    def redact(value):
        if isinstance(value, str):
            return pattern.sub(lambda match: names[match[0]], value)
        if isinstance(value, dict):
            return {redact(key): redact(item) for key, item in value.items()}
        if isinstance(value, list):
            return [redact(item) for item in value]
        if isinstance(value, tuple):
            return tuple(map(redact, value))
        return value

    return redact(report)


def render_html(report: dict, *, hide_folders=None) -> str:
    """Show the assessment as printable HTML, hiding the current and home folders by default."""
    report = _display_report(report, hide_folders)
    assay_id, status, method = _text(report["assay_id"]), report["status"], report["method_report"]
    passed = status == "pass"
    checked = status
    if passed:
        purity = report["purity"]
        percent, _ = _purity_and_standard_text(method)
        expanded = _expanded_text(method["U"])
        if method["U"] is not None:
            expanded += f' with k = {method["k"]["value"]:g}'
        checked = f'{_text(purity["value"])} {purity["unit"]} ({percent}, {expanded})'
    lines = ["<!DOCTYPE html>", '<html lang="en">', "<head>", '<meta charset="utf-8">',
             '<meta name="viewport" content="width=device-width, initial-scale=1">',
             f"<title>{escape(assay_id)} | Spinivo assessment</title>", f"<style>{_STYLE}</style>",
             "</head>", "<body>", "<header>", f"<h1>{escape(assay_id)}</h1>",
             _paragraph(f'Method: {report["method"]} · Mode: {report["mode"]} · Status: {status}'), "</header>",
             "<main>", f'<p class="purity">Checked purity: {escape(checked)}</p>']
    if report["reasons"]:
        lines.extend(_list("Blocking reasons", report["reasons"]))
    if report["not_assessed"]:
        lines.extend(_list(f'Not assessed in {report["mode"]} mode', report["not_assessed"]))
    if method["prerequisites"]:
        lines.append("<h2>Checks</h2>")
        lines.extend(_table("Prerequisites", ("Check", "Status", "Reason"),
                            ((check["name"], check["status"], check["reason"]) for check in method["prerequisites"])))
    if report.get("baseline_shape") is not None:
        lines.extend(("<h2>Baseline shape</h2>", _paragraph(baseline_shape_text(report["baseline_shape"]))))
    if report["runs"]:
        lines.append("<h2>Runs</h2>")
        for label, run in report["runs"].items():
            lines.extend(_run(label, run))
    pairs = report.get("capture_pairs")
    if pairs:
        lines.append("<h2>Capture pairs</h2>")
        lines.extend(_table("Capture pairs", ("Signals", *_CHECK_COLUMNS, "Threshold version"),
                            (("/".join(pair["signals"]), *_check_row(pair), pair["threshold_version"])
                             for pair in pairs)))
    lines.extend(_method_report(method, report["method"] == "external_calibration", passed))
    if report["profiles"]:
        lines.append("<h2>Profiles</h2>")
        for name, profile in report["profiles"].items():
            lines.append(_paragraph(f"{name} profile: {_pairs(profile)}"))
    lines.extend(("</main>", "<footer>", _paragraph(f'{report["schema"]}, spinivo {report["spinivo_version"]}'),
                  _paragraph("Spinivo reports what the evidence supports. It makes no general accuracy claim."),
                  "</footer>", "</body>", "</html>"))
    return "\n".join(lines)
