"""Render assay reports, CRM comparisons and raw FID metadata on the command line."""

import argparse
import os
import sys
from dataclasses import fields
from pathlib import Path

from . import assay_record, assessment_report, crm_comparison, html_report, readers

_INSPECT_SCHEMA = "spinivo.inspect/1"


def _escape_controls(message):
    return "".join(char if char.isprintable() else ascii(char)[1:-1] for char in message)


class _ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        super().error(_escape_controls(message))


def _parser():
    parser = _ArgumentParser(prog="spinivo", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    assess = commands.add_parser("assess", help="assess a JSON assay record")
    assess.add_argument("record", metavar="RECORD", type=Path)
    assess.add_argument("--format", choices=("text", "json", "html"), default="text")
    compare = commands.add_parser("compare", help="compare a checked assay with a CRM certificate")
    compare.add_argument("comparison", metavar="COMPARISON", type=Path)
    inspect = commands.add_parser("inspect", help="show raw FID acquisition metadata")
    inspect.add_argument("path", metavar="PATH", type=Path)
    for command in (compare, inspect):
        command.add_argument("--format", choices=("text", "json"), default="text")
    return parser


def _write_line(stream, message):
    encoding = getattr(stream, "encoding", None)
    if encoding:
        message = message.encode(encoding, errors="backslashreplace").decode(encoding)
    stream.write(message.rstrip("\n") + "\n")


def _diagnostic(message):
    _write_line(sys.stderr, f"spinivo: {_escape_controls(message)}")


def _assess(path, format_name):
    try:
        record = assay_record.load_record(path)
    except (assay_record.RecordError, OSError, RecursionError) as error:
        _diagnostic(f"cannot load record: {error}")
        return 3
    try:
        result = assay_record.assess_record(record)
    except readers.UnsupportedInput as error:
        _diagnostic(f"{error.reason}: {error}")
        return 1
    report = assessment_report.assessment_report(result, record if record.mode == "calculation_only" else None)
    if format_name == "html":
        page = html_report.render_html(report, hide_folders=(path.parent, Path.cwd(), Path.home())) + "\n"
        # Keep the document's bytes consistent with its UTF-8 declaration.
        buffer = getattr(sys.stdout, "buffer", None)
        if buffer is None:
            sys.stdout.write(page)
        else:
            buffer.write(page.encode("utf-8", errors="xmlcharrefreplace"))
        return 0
    render = assessment_report.render_json if format_name == "json" else assessment_report.render_text
    _write_line(sys.stdout, render(report))
    return 0


def _compare(path, format_name):
    try:
        report = crm_comparison.compare(path)
    except (assay_record.RecordError, OSError, RecursionError) as error:
        _diagnostic(f"cannot load record: {error}")
        return 3
    except readers.UnsupportedInput as error:
        _diagnostic(f"{error.reason}: {error}")
        return 1
    render = crm_comparison.render_json if format_name == "json" else crm_comparison.render_text
    _write_line(sys.stdout, render(report))
    return 0


def _inspect_text(report):
    sources = report.get("source_fields", {})
    lines = []
    for field, value in report.items():
        if field in ("schema", "source_fields"):
            continue
        if field == "missing":
            value = ", ".join(value) or "none"
        source = f" ({sources[field]})" if field in sources else ""
        lines.append(f"{field}: {'missing' if value is None else value}{source}")
    return "\n".join(lines)


def _inspect(path, format_name):
    report = {"schema": _INSPECT_SCHEMA}
    try:
        raw = readers.read(path)
    except readers.UnsupportedInput as error:
        report.update(reason=error.reason, detail=error.detail, source=str(error.source))
        _diagnostic(f"{error.reason}: {error}")
    except OSError as error:
        _diagnostic(f"cannot inspect input: {error}")
        return 1
    else:
        report.update((field.name, getattr(raw, field.name)) for field in fields(raw)
                      if field.name not in ("fid", "digital_filter", "source_fields"))
        report["source"] = str(raw.source)
        report["source_fields"] = {field: source for field, source in raw.source_fields.items() if field in report}
    render = assessment_report.render_json if format_name == "json" else _inspect_text
    _write_line(sys.stdout, render(report))
    return 1 if "reason" in report else 0


def _run(argv):
    try:
        args = _parser().parse_args(argv)
    except SystemExit as error:
        return error.code
    if args.command == "assess":
        return _assess(args.record, args.format)
    if args.command == "compare":
        return _compare(args.comparison, args.format)
    return _inspect(args.path, args.format)


def main(argv=None) -> int:
    """Run a Spinivo command and return its exit status."""
    try:
        status = _run(argv)
        sys.stdout.flush()
    except BrokenPipeError:
        # Replace the closed pipe so Python's shutdown flush cannot hit it again.
        with open(os.devnull, "wb") as sink:
            os.dup2(sink.fileno(), sys.stdout.fileno())
        return 1
    return status
