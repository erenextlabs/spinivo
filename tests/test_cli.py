import json
import os
import subprocess
import sys
import unittest
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from spinivo import RecordError, UnsupportedInput, assess_record, load_record, read
from spinivo.assessment_report import assessment_report, render_json, render_text
from spinivo.cli import main
from tests import test_assay_record as records
from tests.test_readers import BRUKER, JEOL, JEOL_ARRAY, PROJECT


def write_synthetic_bruker(path, raw):
    """Store the assay test's synthetic trace so both CLI entry points can read it."""
    path.mkdir()
    parameters = {
        "NUC1": "<1H>", "AQ_mod": 3, "TD": 2 * raw.complex_points, "BYTORDA": 0, "DTYPA": 2,
        "SW_h": raw.spectral_width_hz, "SFO1": raw.observe_frequency_mhz, "BF1": raw.observe_frequency_mhz,
        "O1": raw.carrier_ppm * raw.observe_frequency_mhz, "PULPROG": f"<{raw.pulse_program}>",
        "D": f"(0..1)\n0 {raw.relaxation_delay_s}", "P": f"(0..1)\n0 {raw.pulse_width_us}",
        "GRPDLY": 1e-9, "NC": 0, "DS": 0, "NS": raw.scans or 32,
        "RG": raw.receiver_gain or 30, "TE": raw.temperature_k or 298.15,
    }
    (path / "acqus").write_text("".join(f"##${key}= {value}\n" for key, value in parameters.items()) + "##END=\n",
                               encoding="utf-8")
    (path / "fid").write_bytes(raw.fid.astype("<c16").tobytes())


class CLITests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(TemporaryDirectory()))

    def run_cli(self, *args, env=None, **kwargs):
        return subprocess.run(
            [sys.executable, "-m", "spinivo", *map(str, args)], cwd=PROJECT, timeout=30,
            env={**os.environ, "PYTHONPATH": "src", **(env or {})}, **kwargs,
        )

    def invoke(self, *args):
        argv = list(map(str, args))
        stdout, stderr = StringIO(), StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = main(argv)
        result = code, stdout.getvalue(), stderr.getvalue()
        completed = self.run_cli(*argv, capture_output=True, text=True)
        self.assertEqual((completed.returncode, completed.stdout, completed.stderr), result)
        for output in result[1:]:
            if output:
                self.assertTrue(output.endswith("\n"))
                self.assertFalse(output.endswith("\n\n"))
        return result

    def write_record(self, data):
        path = self.root / "record.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        return path

    def assert_assessment(self, path, status, *, purity_available=False):
        record = load_record(path)
        report = assessment_report(assess_record(record), record)
        self.assertEqual(report["status"], status, report["reasons"])
        self.assertEqual(report["method_report"]["purity"]["value"] is not None, purity_available)
        self.assertEqual(self.invoke("assess", path), (0, render_text(report) + "\n", ""))
        code, stdout, stderr = self.invoke("assess", path, "--format", "json")
        self.assertEqual((code, stdout, stderr), (0, render_json(report) + "\n", ""))
        self.assertEqual(json.loads(stdout, parse_constant=self.fail), report)

    def test_complete_spectral_assessments(self):
        spectral = records.SpectralRecordTests
        spectral.setUpClass()
        for name, raw in spectral.raws.items():
            write_synthetic_bruker(self.root / name, raw)
        for original in (spectral.internal, spectral.external):
            data = deepcopy(original)
            if "pair" in data:
                data["pair"]["source"] = str(records.RECORDS / data["pair"]["source"])
                data["scale"]["storage_scales"] = ["NC=0", "NC=0"]
            with self.subTest(mode=data["mode"]):
                self.assert_assessment(self.write_record(data), "pass", purity_available=True)

    def test_committed_assessments(self):
        for name in (records.JEOL, records.EXTERNAL, records.ENTERED):
            with self.subTest(record=name):
                self.assert_assessment(records.RECORDS / f"{name}.json", "cannot be assessed",
                                       purity_available=name == records.ENTERED)

    def test_assessments_with_missing_evidence(self):
        for name in (records.JEOL, records.EXTERNAL, records.ENTERED):
            with self.subTest(record=name):
                data = records.record_data(name, ("N_a",), None)
                for role in ("run", "sample", "calibrant"):
                    run = data.get(role, {})
                    if "file" in run:
                        run["file"] = str(records.RECORDS / run["file"])
                if "pair" in data:
                    data["pair"]["source"] = str(records.RECORDS / data["pair"]["source"])
                self.assert_assessment(self.write_record(data), "cannot be assessed")

    def test_assess_unsupported_spectrum_exits_one(self):
        with self.assertRaises(UnsupportedInput) as caught:
            read(JEOL_ARRAY)
        error = caught.exception
        for name, role in ((records.JEOL, "run"), (records.EXTERNAL, "sample"), (records.EXTERNAL, "calibrant")):
            with self.subTest(role=role):
                data = records.record_data(name)
                for run in ("run", "sample", "calibrant"):
                    if run in data:
                        data[run]["file"] = str(JEOL_ARRAY if run == role else records.RECORDS / data[run]["file"])
                self.assertEqual(self.invoke("assess", self.write_record(data)),
                                 (1, "", f"spinivo: {error.reason}: {error}\n"))

    def test_assess_reader_error_still_renders_withheld_report(self):
        path = records.RECORDS / f"{records.JEOL}.json"
        stdout, stderr = StringIO(), StringIO()
        with (patch("spinivo.pipeline.read", side_effect=OSError("corrupt data")),
              redirect_stdout(stdout), redirect_stderr(stderr)):
            code = main(["assess", str(path), "--format=json"])
        self.assertEqual((code, stderr.getvalue()), (0, ""))
        report = json.loads(stdout.getvalue())
        self.assertEqual(report["status"], "cannot be assessed")
        self.assertIn("reader: corrupt data", report["reasons"])

    def test_malformed_assay_records(self):
        for name in (records.JEOL, records.EXTERNAL, records.ENTERED):
            with self.subTest(record=name):
                path = self.write_record(records.record_data(name, ("N_a",), "two"))
                with self.assertRaises(RecordError) as caught:
                    load_record(path)
                self.assertEqual(self.invoke("assess", path, "--format=json"),
                                 (3, "", f"spinivo: cannot load record: {caught.exception}\n"))

    def test_missing_run_objects_are_malformed(self):
        for name, role in ((records.JEOL, "run"), (records.EXTERNAL, "sample"), (records.EXTERNAL, "calibrant")):
            with self.subTest(role=role):
                data = records.record_data(name)
                del data[role]
                path = self.write_record(data)
                self.assertEqual(self.invoke("assess", path, "--format=json"),
                                 (3, "", f"spinivo: cannot load record: {role}.file: must name the data file\n"))

    def test_stdout_escapes_unencodable_characters(self):
        record = self.write_record(records.record_data(records.ENTERED, ("assay_id",), "\ud800"))
        fid = self.root / "café.jdf"
        fid.write_bytes(JEOL.read_bytes())
        for command, path, field, value, escaped in (
            ("assess", record, "assay_id", "\ud800", r"\ud800"),
            ("inspect", fid, "source", str(fid), r"caf\xe9.jdf"),
        ):
            for format_name in ("text", "json"):
                with self.subTest(command=command, format=format_name):
                    result = self.run_cli(command, path, "--format", format_name, capture_output=True,
                                          encoding="ascii", env={"PYTHONIOENCODING": "ascii:strict"})
                    self.assertEqual((result.returncode, result.stderr), (0, ""))
                    self.assertTrue(result.stdout.endswith("\n"))
                    if format_name == "json":
                        self.assertEqual(json.loads(result.stdout)[field], value)
                    else:
                        self.assertIn(escaped, result.stdout)

    def test_stderr_escapes_control_characters(self):
        key = "bad\nkey\r\x1b[31m\u0085\u202e"
        escaped = r"bad\nkey\r\x1b[31m\x85\u202e"
        path = self.write_record({**records.record_data(records.ENTERED), key: None})
        self.assertEqual(self.invoke("assess", path),
                         (3, "", f"spinivo: cannot load record: {escaped}: unknown field\n"))
        code, stdout, stderr = self.invoke("inspect", self.root / key, "--format=json")
        self.assertEqual(code, 1)
        self.assertEqual(json.loads(stdout)["source"], str(self.root / key))
        self.assertIn(escaped, stderr)
        self.assertEqual(len(stderr.splitlines()), 1)
        code, stdout, stderr = self.invoke("assess", path, "--" + key)
        self.assertEqual((code, stdout), (2, ""))
        self.assertTrue(stderr.endswith(f"error: unrecognized arguments: --{escaped}\n"))

    def test_closed_stdout_pipe(self):
        for buffered in (True, False):
            with self.subTest(buffered=buffered):
                read_fd, write_fd = os.pipe()
                os.close(read_fd)
                with os.fdopen(write_fd, "wb") as sink:
                    result = self.run_cli("inspect", BRUKER, "--format=json", stdout=sink, stderr=subprocess.PIPE,
                                          env={"PYTHONUNBUFFERED": "0" if buffered else "1"})
                self.assertEqual((result.returncode, result.stderr), (1, b""))

    def test_record_read_errors(self):
        invalid = self.root / "invalid.json"
        for content in (b"{", b"\xff", b"[" * 2000):
            invalid.write_bytes(content)
            with self.subTest(content=content[:8]):
                code, stdout, stderr = self.invoke("assess", invalid)
                self.assertEqual((code, stdout), (3, ""))
                self.assertTrue(stderr.startswith("spinivo: cannot load record: "))
                self.assertEqual(len(stderr.splitlines()), 1)
        for path in (self.root / "missing.json", self.root):
            with self.subTest(path=path):
                with self.assertRaises(OSError) as caught:
                    load_record(path)
                self.assertEqual(self.invoke("assess", path),
                                 (3, "", f"spinivo: cannot load record: {caught.exception}\n"))

    def test_inspect_committed_fids(self):
        for path in (JEOL, BRUKER):
            with self.subTest(path=path):
                raw = read(path)
                expected = {key: value for key, value in vars(raw).items()
                            if value is None or isinstance(value, (int, float, str))}
                expected.pop("digital_filter", None)
                expected["source"] = str(raw.source)
                sources = {key: value for key, value in raw.source_fields.items() if key in expected}
                code, stdout, stderr = self.invoke("inspect", path, "--format", "json")
                self.assertEqual((code, stderr), (0, ""))
                self.assertEqual(json.loads(stdout, parse_constant=self.fail),
                                 {"schema": "spinivo.inspect/1", **expected, "source_fields": sources,
                                  "missing": list(raw.missing)})
                code, stdout, stderr = self.invoke("inspect", path, "--format", "text")
                self.assertEqual((code, stderr), (0, ""))
                expected["missing"] = ", ".join(raw.missing) or "none"
                self.assertCountEqual(stdout.splitlines(), [
                    f"{key}: {'missing' if value is None else value}"
                    + (f" ({sources[key]})" if key in sources else "") for key, value in expected.items()
                ])

    def test_inspect_unsupported_inputs(self):
        truncated = self.root / "truncated.jdf"
        truncated.write_bytes(JEOL.read_bytes()[:128])
        missing_acqus = self.root / "missing acqus"
        missing_acqus.mkdir()
        (missing_acqus / "fid").touch()
        bad_header = self.root / "Bruker: incomplete"
        bad_header.mkdir()
        (bad_header / "fid").touch()
        (bad_header / "acqus").touch()
        for path, reason, source in ((truncated, "truncated", truncated), (JEOL_ARRAY, "arrayed_or_2d", JEOL_ARRAY),
                                     (missing_acqus, "missing_acqus", missing_acqus),
                                     (bad_header, "truncated", bad_header / "acqus")):
            with self.subTest(path=path):
                with self.assertRaises(UnsupportedInput) as caught:
                    read(path)
                self.assertEqual(caught.exception.reason, reason)
                message = str(caught.exception)
                self.assertTrue(message.startswith(f"{source}: "))
                detail = message.removeprefix(f"{source}: ")
                self.assertEqual((caught.exception.source, caught.exception.detail), (source, detail))
                stderr = f"spinivo: {reason}: {message}\n"
                code, stdout, actual_stderr = self.invoke("inspect", path, "--format", "json")
                self.assertEqual((code, actual_stderr), (1, stderr))
                self.assertEqual(json.loads(stdout, parse_constant=self.fail), {
                    "schema": "spinivo.inspect/1", "reason": reason, "detail": detail, "source": str(source),
                })
                self.assertEqual(self.invoke("inspect", path),
                                 (1, f"reason: {reason}\ndetail: {detail}\nsource: {source}\n", stderr))

    def test_usage_and_help(self):
        for args in ((), ("assess",), ("unknown",), ("inspect", JEOL, "--format", "xml")):
            with self.subTest(args=args):
                code, stdout, stderr = self.invoke(*args)
                self.assertEqual((code, stdout), (2, ""))
                self.assertIn("usage: spinivo", stderr)
                self.assertIn("error:", stderr)
        code, stdout, stderr = self.invoke("--help")
        self.assertEqual((code, stderr), (0, ""))
        self.assertIn("usage: spinivo", stdout)
