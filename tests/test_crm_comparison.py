import json
import unittest
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from io import StringIO
from math import sqrt
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from spinivo import RecordError, __version__, assess_record, load_record
from spinivo.assessment_report import assessment_report
from spinivo.cli import main
from spinivo.crm_comparison import compare, render_json, render_text
from tests import test_assay_record as records, test_cli

AGREEMENT, LIMIT = "uncertainty_agreement", "lab_limit"
WITHHELD = "cannot be assessed"


class CRMComparisonTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(TemporaryDirectory()))
        assay_folder = self.root / "assays"
        assay_folder.mkdir()
        spectral = records.SpectralRecordTests
        spectral.setUpClass()
        test_cli.write_synthetic_bruker(assay_folder / "internal.jdf", spectral.raws["internal.jdf"])
        self.assay_path = assay_folder / "assay.json"
        self.assay_data = deepcopy(spectral.internal)
        self.assay_path.write_text(json.dumps(self.assay_data), encoding="utf-8")
        self.path = self.root / "comparison.json"
        self.data = {
            "schema": "spinivo.crm-comparison/1", "comparison_id": "C01", "assay_record": "assays/assay.json",
            "crm": {"id": "CRM01", "lot": "L01", "certificate_source": "Certificate CRM01 L01",
                    "mass_fraction_pct": 99.5, "standard_u_pp": 0.5},
            "independence": {"independent_of_calibration": True, "basis": "Separate reference and calibration inputs"},
            "comparison_k": {"value": 2, "basis": "Lab uncertainty protocol"},
            "lab_limit_pp": {"value": 1.0, "basis": "Lab protocol fixed before measurement"},
        }

    run_cli = test_cli.CLITests.run_cli
    invoke = test_cli.CLITests.invoke

    def report(self):
        self.path.write_text(json.dumps(self.data), encoding="utf-8")
        return compare(self.path)

    def assay_report(self):
        record = load_record(self.assay_path)
        return assessment_report(assess_record(record), record)

    def assert_statuses(self, report, agreement, limit):
        self.assertEqual((report[AGREEMENT]["status"], report[LIMIT]["status"]), (agreement, limit))

    def test_hand_worked_and_expanded_uncertainty(self):
        assay = self.assay_report()
        self.assertEqual(assay["status"], "pass", assay["reasons"])
        measured = assay["purity"]["value"] * 100
        measurement_u = assay["method_report"]["u_c"] * 100
        error = measured - 99.5
        error_u = sqrt(measurement_u**2 + 0.5**2)
        standard = self.report()
        self.assert_statuses(standard, "pass", "pass")
        self.assertEqual(standard["schema"], "spinivo.crm-comparison-report/1")
        self.assertEqual(standard["spinivo_version"], __version__)
        self.assertEqual(standard["comparison_id"], "C01")
        self.assertEqual(standard["crm"], {"id": "CRM01", "lot": "L01", "certificate_source": "Certificate CRM01 L01"})
        self.assertEqual(standard["assay"], {key: assay[key] for key in ("assay_id", "method", "status", "reasons")})
        for field, expected in (
            ("measured_purity_pct", measured), ("reference_purity_pct", 99.5),
            ("measurement_standard_u_pp", measurement_u), ("reference_standard_u_pp", 0.5),
            ("error_pp", error), ("absolute_error_pp", abs(error)), ("error_standard_u_pp", error_u),
            ("k", 2), ("comparison_u_pp", 2 * error_u), ("standard_u_ratio", abs(error) / (2 * error_u)),
        ):
            with self.subTest(field=field):
                self.assertAlmostEqual(standard[field], expected, places=12)
        self.data["crm"].update(standard_u_pp=None, expanded_u_pp=1.5, coverage_factor=3)
        self.assertEqual(self.report(), standard)

    def test_each_verdict_boundary(self):
        assay = self.assay_report()
        # Both thresholds are 0.4 pp, which subtraction at 99 % can round slightly above.
        assay["method_report"]["u_c"] = 0
        assay["purity"]["value"] = 0.99
        self.data["crm"]["standard_u_pp"] = 0.4
        self.data["comparison_k"]["value"] = 1
        self.data["lab_limit_pp"]["value"] = 0.4
        with patch("spinivo.crm_comparison.assessment_report.assessment_report", return_value=assay):
            for sign in (-1, 1):
                for absolute, expected in ((0.39, "pass"), (0.4, "pass"), (0.41, "fail")):
                    with self.subTest(sign=sign, absolute=absolute):
                        self.data["crm"]["mass_fraction_pct"] = 99 - sign * absolute
                        report = self.report()
                        self.assertAlmostEqual(report["error_pp"], sign * absolute)
                        self.assert_statuses(report, expected, expected)
            self.data["lab_limit_pp"]["value"] = 1
            self.assert_statuses(self.report(), "fail", "pass")
            self.data["comparison_k"]["value"] = 2
            self.data["lab_limit_pp"]["value"] = 0.4
            self.assert_statuses(self.report(), "pass", "fail")

    def test_remove_and_restore_evidence(self):
        original = self.report()
        for field, reason, withheld in (
            ("assay_record", "Missing assay_record", (AGREEMENT, LIMIT)),
            ("crm.mass_fraction_pct", "Missing crm.mass_fraction_pct", (AGREEMENT, LIMIT)),
            ("independence.independent_of_calibration", "must be declared true", (AGREEMENT, LIMIT)),
            ("independence.basis", "Missing independence.basis", (AGREEMENT, LIMIT)),
            ("comparison_k.value", "Missing comparison_k.value", (AGREEMENT,)),
            ("comparison_k.basis", "Missing comparison_k.basis", (AGREEMENT,)),
            ("crm.standard_u_pp", "Missing crm.standard_u_pp", (AGREEMENT,)),
            ("lab_limit_pp.value", "Missing lab_limit_pp.value", (LIMIT,)),
            ("lab_limit_pp.basis", "Missing lab_limit_pp.basis", (LIMIT,)),
        ):
            with self.subTest(field=field):
                parts = field.split(".")
                parent = self.data[parts[0]] if len(parts) == 2 else self.data
                value = parent.pop(parts[-1])
                self.assert_gap(self.report(), original, reason, withheld)
                parent[parts[-1]] = value
                self.assertEqual(self.report(), original)

        crm = self.data["crm"]
        crm.update(standard_u_pp=None, expanded_u_pp=1.5, coverage_factor=3)
        for field in ("expanded_u_pp", "coverage_factor"):
            with self.subTest(field=field):
                old, crm[field] = crm[field], None
                self.assert_gap(self.report(), original, f"crm.{field}", (AGREEMENT,))
                crm[field] = old
                self.assertEqual(self.report(), original)

        assay = self.assay_report()
        assay["method_report"]["u_c"] = None
        with patch("spinivo.crm_comparison.assessment_report.assessment_report", return_value=assay):
            self.assert_gap(self.report(), original, "Missing measurement standard uncertainty", (AGREEMENT,))
        self.assertEqual(self.report(), original)

    def test_explicit_non_independence_withholds_both_verdicts(self):
        original = self.report()
        independence = self.data["independence"]
        independence["independent_of_calibration"] = False
        self.assert_gap(self.report(), original, "independence.independent_of_calibration must be declared true",
                        (AGREEMENT, LIMIT))
        independence["independent_of_calibration"] = True
        self.assertEqual(self.report(), original)

    def assert_gap(self, report, original, reason, withheld):
        for verdict in (AGREEMENT, LIMIT):
            if verdict in withheld:
                self.assertEqual(report[verdict]["status"], WITHHELD)
                self.assertIn(reason, report[verdict]["reason"])
            else:
                self.assertEqual(report[verdict], original[verdict])

    def test_cli_matches_library_and_retains_assay_reasons(self):
        passing = self.report()
        for withheld in (False, True):
            with self.subTest(withheld=withheld):
                if withheld:
                    self.assay_data["N_a"] = None
                    self.assay_path.write_text(json.dumps(self.assay_data), encoding="utf-8")
                report = self.report()
                if withheld:
                    assay = self.assay_report()
                    self.assertTrue(assay["reasons"])
                    self.assertEqual(report["assay"], {key: assay[key] for key in passing["assay"]})
                    self.assertIsNone(report["measured_purity_pct"])
                    self.assertIsNone(report["error_pp"])
                    self.assert_gap(report, passing, "A checked pass is required", (AGREEMENT, LIMIT))
                self.assertEqual(self.invoke("compare", self.path), (0, render_text(report) + "\n", ""))
                self.assertEqual(self.invoke("compare", self.path, "--format", "json"),
                                 (0, render_json(report) + "\n", ""))
                self.assertEqual(json.loads(render_json(report), parse_constant=self.fail), report)
                self.assertIn("Neither verdict qualifies a method or releases a batch.", render_text(report))

    def test_text_escapes_certificate_source_and_assay_reason(self):
        source = "Certificate\nrevision\t2\x1b[31m"
        self.data["crm"]["certificate_source"] = source
        self.assay_data["signals"]["analyte"][0].update(name="Peak\nname\x1b[31m", bounds_ppm=None)
        self.assay_path.write_text(json.dumps(self.assay_data), encoding="utf-8")
        report = self.report()
        rendered = render_text(report)
        lines = rendered.splitlines()
        self.assertIn(r"Certificate source: Certificate\nrevision\t2\x1b[31m", lines)
        self.assertIn(r"Assay reason: Peak\nname\x1b[31m analyte bounds: Declare two distinct finite ppm bounds", lines)
        self.assertNotIn("\x1b", rendered)
        payload = json.loads(render_json(report))
        self.assertEqual(payload["crm"]["certificate_source"], source)
        self.assertIn("Peak\nname\x1b[31m analyte bounds: Declare two distinct finite ppm bounds",
                      payload["assay"]["reasons"])

    def test_unsupported_reader_exits_one(self):
        self.assay_data["run"]["file"] = str(test_cli.JEOL_ARRAY)
        self.assay_path.write_text(json.dumps(self.assay_data), encoding="utf-8")
        self.path.write_text(json.dumps(self.data), encoding="utf-8")
        result = self.invoke("compare", self.path)
        self.assertEqual(result, self.invoke("assess", self.assay_path))
        self.assertEqual(result[0], 1)

    def test_malformed_inputs_exit_three(self):
        for field, value in (
            ("schema", "spinivo.assay/1"), ("comparison_id", 1), ("assay_record", []),
            ("assay_record", "/assay.json"), ("crm", []), ("crm.lot", False),
            ("crm.mass_fraction_pct", True), ("crm.mass_fraction_pct", 0), ("crm.mass_fraction_pct", 101),
            ("crm.standard_u_pp", -0.1), ("crm.standard_u_pp", float("inf")), ("crm.standard_u_pp", 10**400),
            ("crm.expanded_u_pp", 1), ("crm.coverage_factor", 0), ("comparison_k.value", 0),
            ("independence.independent_of_calibration", 1), ("independence.basis", {}),
            ("lab_limit_pp.value", -1), ("lab_limit_pp.basis", []), ("comparison_k.value", float("nan")),
        ):
            with self.subTest(field=field, value=value):
                data = deepcopy(self.data)
                parts = field.split(".")
                parent = data[parts[0]] if len(parts) == 2 else data
                parent[parts[-1]] = value
                self.path.write_text(json.dumps(data), encoding="utf-8")
                self.assert_bad_record("crm:" if field == "crm.expanded_u_pp" else field + ":")
        for content in ("{", "null", '{"schema": "spinivo.crm-comparison/1", "comparison_k": {"value": 1e-400}}'):
            with self.subTest(content=content):
                self.path.write_text(content, encoding="utf-8")
                self.assert_bad_record("")

    def assert_bad_record(self, field):
        with self.assertRaises(RecordError) as caught:
            compare(self.path)
        self.assertTrue(str(caught.exception).startswith(field))
        stdout, stderr = StringIO(), StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = main(["compare", str(self.path)])
        self.assertEqual((code, stdout.getvalue(), stderr.getvalue()),
                         (3, "", f"spinivo: cannot load record: {caught.exception}\n"))
