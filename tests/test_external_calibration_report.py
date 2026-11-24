import json
import unittest
from dataclasses import asdict, replace
from math import inf, isclose, nan
from pathlib import Path

from spinivo.external_calibration import ExternalRun, calculate_external_calibration, link_external_checks
from spinivo.external_calibration_report import purity_report, render_json, render_text
from spinivo.internal_standard import CALCULATION_ONLY, Quantity
from spinivo.internal_standard_report import render_text as render_internal_text
from tests.integral_fixtures import load_cases, make_raw
from tests.test_external_calibration import PREPARATION, QUANTITIES, load_assay_case
from tests.test_integral_report import fixture_result
from tests.test_internal_standard import acquisition_checks


class ExternalCalibrationReportTests(unittest.TestCase):
    def setUp(self):
        self.assay, self.expected = load_assay_case()

    def report(self, assay):
        report = purity_report(calculate_external_calibration(assay), assay)
        self.assertEqual(json.loads(json.dumps(report, allow_nan=False)), report)
        return json.loads(render_json(report))

    def assertRelative(self, actual, expected):
        self.assertTrue(isclose(actual, expected, rel_tol=1e-9), (actual, expected))

    def test_separate_external_result_cannot_report_pass(self):
        result = replace(calculate_external_calibration(self.assay), status="pass", reasons=())
        with self.assertRaisesRegex(ValueError, "checked pass report"):
            purity_report(result, self.assay)

    def test_missing_run_sources_are_null(self):
        report = self.report(self.assay)
        for role in ("sample", "calibrant"):
            self.assertIsNone(report[role]["header"])
            self.assertIsNone(report[role]["entered"])

    def test_linked_run_without_entered_values_is_null(self):
        blank = ExternalRun(*[None] * 12)
        assay = replace(self.assay, _linked_runs=(("sample", blank, self.assay.sample),
                                                  ("calibrant", replace(blank, scans=32), self.assay.calibrant)))
        report = self.report(assay)
        self.assertEqual((report["sample"]["entered"], report["calibrant"]["entered"]), (None, {"scans": 32}))

    def test_report_uses_calculation_record_and_keeps_known_unscaled_areas(self):
        result = calculate_external_calibration(self.assay)
        incomplete = replace(self.assay, pair=None, S=None, sample=replace(self.assay.sample, reader_scale=None))
        self.assertEqual(purity_report(result, incomplete)["pair"], asdict(self.assay.pair))
        self.assertEqual(purity_report(result, incomplete), purity_report(result))
        report = self.report(incomplete)
        self.assertIsNone(report["purity"]["value"])
        self.assertIsNone(report["J_a"])
        run = report["sample"]
        self.assertEqual((run["area"], run["scans"], run["per_scan_area"]), (0.55648986, 32, 0.55648986 / 32))
        self.assertIn("sample area = 0.55648986, scans = 32, area per scan = 0.017390308", render_text(report))

    def test_saved_version_one_and_two_reports_show_only_their_fields(self):
        # Serialized output from the archived writers for the same certified assay.
        v1 = json.loads("""{"J_a": 1.0, "J_s": 1.0, "K_covers": [], "R": 1.0, "U": 0.02, "U_reason": null, "assay_id": "Simple",
            "budget": [{"contribution": 0.0, "name": "R", "sensitivity": 1.0, "share": 0.0}, {"contribution": 0.0, "name": "c_s", "sensitivity": 1.0, "share": 0.0}, {"contribution": 0.0, "name": "M_a", "sensitivity": 0.01, "share": 0.0}, {"contribution": 0.0, "name": "V_a", "sensitivity": 0.1, "share": 0.0}, {"contribution": -0.0, "name": "m_a", "sensitivity": -1.0, "share": 0.0}, {"contribution": 0.01, "name": "K", "sensitivity": 1.0, "share": 1.0}, {"contribution": 0.0, "name": "D", "sensitivity": 1.0, "share": 0.0}], "c_a": 1.0, "c_s": 1,
            "calibrant": {"area": 32, "area_mode": "summed", "decoupling": "off", "preparation_id": "s", "probe": "P", "pulse_sequence": "zg", "receiver_gain": 1, "run_id": "s", "scans": 32, "solvent": "DMSO", "temperature": 25, "tube": "5 mm"}, "concentration_unit": "mmol/L", "covariances": [],
            "inputs": {"D": {"source": "drift log", "u": 0, "unit": "ratio", "value": 1}, "K": {"source": "response log", "u": 0.01, "unit": "ratio", "value": 1}, "M_a": {"source": "formula", "u": 0, "unit": "g/mol", "value": 100}, "N_a": {"source": null, "u": null, "unit": "count", "value": 1}, "N_s": {"source": null, "u": null, "unit": "count", "value": 1}, "R": {"source": "repeatability", "u": 0.0, "unit": "ratio", "value": 1.0}, "R_relative_u": {"source": "repeatability", "u": null, "unit": "fraction", "value": 0}, "V_a": {"source": "volume", "u": 0, "unit": "mL", "value": 10}, "c_s": {"source": "certificate", "u": 0, "unit": "mmol/L", "value": 1}, "m_a": {"source": "mass", "u": 0, "unit": "mg", "value": 1}}, "k": {"source": "coverage", "value": 2}, "method": "external calibration", "prerequisites": [], "purity": {"unit": "fraction", "value": 1.0}, "reasons": [],
            "sample": {"area": 32, "area_mode": "summed", "decoupling": "off", "preparation_id": "a", "probe": "P", "pulse_sequence": "zg", "receiver_gain": 1, "run_id": "a", "scans": 32, "solvent": "DMSO", "temperature": 25, "tube": "5 mm"},
            "schema": "spinivo.external_calibration_purity/1", "status": "pass", "u_c": 0.01, "u_c/w": 0.01, "variance": 0.0001}""")
        v2 = json.loads("""{"J_a": 1.0, "J_s": 1.0, "K_covers": [], "R": 1.0, "U": 0.02, "U_reason": null, "assay_id": "Simple",
            "budget": [{"contribution": 0.0, "name": "R", "sensitivity": 1.0, "share": 0.0}, {"contribution": 0.0, "name": "c_s", "sensitivity": 1.0, "share": 0.0}, {"contribution": 0.0, "name": "M_a", "sensitivity": 0.01, "share": 0.0}, {"contribution": 0.0, "name": "V_a", "sensitivity": 0.1, "share": 0.0}, {"contribution": -0.0, "name": "m_a", "sensitivity": -1.0, "share": 0.0}, {"contribution": 0.01, "name": "K", "sensitivity": 1.0, "share": 1.0}, {"contribution": 0.0, "name": "D", "sensitivity": 1.0, "share": 0.0}], "c_a": 1.0, "c_s": 1, "c_s_route": "certified",
            "calibrant": {"area": 32, "area_mode": "summed", "decoupling": "off", "preparation_id": "s", "probe": "P", "pulse_sequence": "zg", "receiver_gain": 1, "run_id": "s", "scans": 32, "solvent": "DMSO", "temperature": 25, "tube": "5 mm"}, "concentration_unit": "mmol/L", "covariances": [],
            "inputs": {"D": {"source": "drift log", "u": 0, "unit": "ratio", "value": 1}, "K": {"source": "response log", "u": 0.01, "unit": "ratio", "value": 1}, "M_a": {"source": "formula", "u": 0, "unit": "g/mol", "value": 100}, "N_a": {"source": null, "u": null, "unit": "count", "value": 1}, "N_s": {"source": null, "u": null, "unit": "count", "value": 1}, "R": {"source": "repeatability", "u": 0.0, "unit": "ratio", "value": 1.0}, "R_relative_u": {"source": "repeatability", "u": null, "unit": "fraction", "value": 0}, "V_a": {"source": "volume", "u": 0, "unit": "mL", "value": 10}, "c_s": {"source": "certificate", "u": 0, "unit": "mmol/L", "value": 1}, "m_a": {"source": "mass", "u": 0, "unit": "mg", "value": 1}}, "k": {"source": "coverage", "value": 2}, "method": "external calibration", "prerequisites": [], "purity": {"unit": "fraction", "value": 1.0}, "reasons": ["calculation only: no spectral evidence"],
            "sample": {"area": 32, "area_mode": "summed", "decoupling": "off", "preparation_id": "a", "probe": "P", "pulse_sequence": "zg", "receiver_gain": 1, "run_id": "a", "scans": 32, "solvent": "DMSO", "temperature": 25, "tube": "5 mm"},
            "schema": "spinivo.external_calibration_purity/2", "status": "cannot be assessed", "u_c": 0.01, "u_c/w": 0.01, "variance": 0.0001}""")
        for report, route in ((v1, ""), (v2, " (certified)")):
            with self.subTest(schema=report["schema"]):
                rendered = render_text(report)
                self.assertIn("sample area = 32, scans = 32\n", rendered)
                self.assertIn("Purity = 100.0 %, u_c = 1.0 %", rendered)
                self.assertIn(f"c_s = 1 mmol/L{route}\n", rendered)
                self.assertIn("  K = 1 (response log)", rendered)
                self.assertNotIn("area per scan", rendered)
                self.assertNotIn("  S =", rendered)

    def test_hand_worked_json_and_text(self):
        report = self.report(self.assay)
        self.assertEqual(report["schema"], "spinivo.external_calibration_purity/3")
        self.assertEqual(report["method"], "external calibration")
        self.assertEqual(report["assay_id"], self.assay.assay_id)
        self.assertEqual(report["status"], "cannot be assessed")
        self.assertEqual(report["purity"]["unit"], "fraction")
        self.assertRelative(report["purity"]["value"], self.expected["purity"])
        self.assertEqual(report["concentration_unit"], "mmol/L")
        self.assertEqual(report["c_s_route"], "prepared")
        for role, area in (("sample", "J_a"), ("calibrant", "J_s")):
            run = getattr(self.assay, role)
            self.assertEqual(report[role], {**asdict(run), "per_scan_area": run.area / run.scans,
                                            "header": None, "entered": None})
            self.assertEqual(report[area], run.area / run.scans)
        self.assertEqual((report["pair"], report["S_gains"]), (asdict(self.assay.pair), [30, 30]))
        self.assertEqual(report["tolerances"], {"temperature": 0.1, "observe_frequency_mhz": 1e-4})
        for name, expected in (("R", "R"), ("c_s", "c_s"), ("c_a", "c_a"), ("variance", "variance"),
                               ("u_c", "standard_u"), ("u_c/w", "relative_standard_u"), ("U", "expanded_u")):
            self.assertRelative(report[name], self.expected[expected])
        self.assertEqual(report["k"], asdict(self.assay.k))
        self.assertIsNone(report["U_reason"])
        self.assertEqual((report["K_covers"], report["prerequisites"], report["reasons"]), ([], [], [CALCULATION_ONLY]))

        inputs = report["inputs"]
        self.assertEqual(set(inputs), {"R", "R_relative_u", "N_a", "N_s", *QUANTITIES, "S"})
        relative_u = self.assay.R_relative_u
        source = self.assay.R_u_source
        self.assertEqual(inputs["R"], {"value": report["R"], "u": report["R"] * relative_u,
                                       "unit": "ratio", "source": source})
        self.assertEqual(inputs["R_relative_u"], {"value": relative_u, "u": None, "unit": "fraction", "source": source})
        for name in ("N_a", "N_s"):
            self.assertEqual(inputs[name], {"value": getattr(self.assay, name), "u": None,
                                            "unit": "count", "source": None})
        for name in (*QUANTITIES, "S"):
            quantity = getattr(self.assay, name)
            self.assertEqual(inputs[name], {"value": quantity.value, "u": quantity.standard_u,
                                            "unit": quantity.unit, "source": quantity.source,
                                            **({"source_ids": list(quantity.source_ids)} if quantity.source_ids else {})})
        self.assertRelative(inputs["K"]["value"], self.expected["K"])

        self.assertEqual({row["name"] for row in report["budget"]}, set(self.expected["sensitivities"]))
        for row in report["budget"]:
            name = row["name"]
            contribution = self.expected["sensitivities"][name] * inputs[name]["u"]
            self.assertRelative(row["sensitivity"], self.expected["sensitivities"][name])
            self.assertRelative(row["contribution"], contribution)
            self.assertRelative(row["share"], contribution**2 / self.expected["variance"])
        self.assertEqual(len(report["covariances"]), len(self.assay.correlations))
        for row, correlation, term in zip(report["covariances"], self.assay.correlations,
                                          self.expected["covariance_terms"]):
            self.assertEqual((row["inputs"], row["r"], row["source"]),
                             (list(correlation.inputs), correlation.r, correlation.source))
            self.assertRelative(row["term"], term)
            self.assertRelative(row["share"], term / self.expected["variance"])
        rendered = render_text(report)
        for line in (CALCULATION_ONLY, "Purity = 97.35 %, u_c = 0.66 %", "U = 1.3 % with k = 2",
                     "c_s = 2.7133857 mmol/L", "c_a = 8.4278497 mmol/L", f"K = 0.99869024 ({self.assay.K.source})",
                     f"S = 0.98 ({self.assay.S.source})", "Contribution", "Variance term", "m_s, m_a"):
            self.assertIn(line, rendered)

    def test_blocked_transfer_keeps_each_run_visible(self):
        report = self.report(replace(self.assay, pair=None, S=None))
        self.assertEqual((report["purity"]["value"], report["pair"], report["inputs"]["S"]["value"]),
                         (None, None, None))
        rendered = render_text(report)
        for line in ("pair record is missing", "S is missing", "equal nominal receiver gains are not scale evidence",
                     "J_a = 0.017390308", "J_s = 0.010959454", "R = 1.586786", "c_s = 2.7133857 mmol/L (prepared)"):
            self.assertIn(line, rendered)

    def test_coverage_rounding_matches_internal_standard(self):
        assay = replace(self.assay, k=replace(self.assay.k, value=2.1))
        report = self.report(assay)
        internal = {**report, "I_a": report["J_a"], "I_s": report["J_s"]}
        expected = next(line for line in render_internal_text(internal).splitlines() if line.startswith("  U = "))
        self.assertEqual(expected, "  U = 1.4 % with k = 2.1")
        self.assertIn(expected, render_text(report).splitlines())

    def test_certified_concentration_without_coverage_factor(self):
        c_s = Quantity(self.expected["c_s"], 0.01, "mmol/L", "Illustrative certificate")
        sample = replace(self.assay.sample, solvent="other solvent")
        assay = replace(self.assay, c_s_route="certified", c_s=c_s, **dict.fromkeys(PREPARATION), sample=sample,
                        K_covers=("solvent",), correlations=(), k=None)
        report = self.report(assay)
        self.assertEqual((report["status"], report["reasons"]), ("cannot be assessed", [CALCULATION_ONLY]))
        self.assertEqual(report["K_covers"], ["solvent"])
        self.assertEqual(report["c_s"], c_s.value)
        self.assertGreater(report["u_c"], 0)
        self.assertEqual(report["inputs"]["c_s"], {"value": c_s.value, "u": c_s.standard_u,
                                                   "unit": c_s.unit, "source": c_s.source})
        self.assertTrue(set(PREPARATION).isdisjoint(report["inputs"]))
        self.assertIsNone(report["k"])
        self.assertIsNone(report["U"])
        self.assertEqual(report["U_reason"], "coverage factor k was not supplied")
        self.assertIn("U unavailable (coverage factor k was not supplied)", render_text(report))

    def test_conflicting_concentration_paths_keep_all_inputs(self):
        c_s = Quantity(self.expected["c_s"], 0.01, "mmol/L", "Illustrative certificate")
        assay = replace(self.assay, c_s=c_s)
        report = self.report(assay)
        self.assertEqual(report["status"], "cannot be assessed")
        self.assertIsNone(report["purity"]["value"])
        self.assertIn("only the declared prepared c_s inputs may be supplied", report["reasons"])
        self.assertIsNone(report["c_s"])
        for name in ("c_s", *PREPARATION):
            quantity = getattr(assay, name)
            self.assertEqual(report["inputs"][name], {"value": quantity.value, "u": quantity.standard_u,
                                                      "unit": quantity.unit, "source": quantity.source,
                                                      **({"source_ids": list(quantity.source_ids)} if quantity.source_ids else {})})

    def test_flagged_linked_region_keeps_per_scan_areas(self):
        original = calculate_external_calibration(self.assay)
        case = next(case for case in load_cases()["cases"] if case["id"] == "clean_singlet")
        raw = make_raw(case)
        clean = fixture_result("clean_singlet")
        inputs = {}
        for role in ("sample", "calibrant"):
            run = getattr(self.assay, role)
            signals = getattr(self.assay, f"{role}_signals")
            region = replace(clean, nominal_area=run.area / len(signals), run=run.run_id)
            if role == "calibrant":
                region = replace(region, status="flag", flags=(("baseline", "distorted"),))
            inputs[f"{role}_raw"] = replace(raw, source=Path(run.run_id), scans=run.scans,
                                            receiver_gain=run.receiver_gain, pulse_program=run.pulse_sequence,
                                            temperature_k=298.15, intensity_scale=1.0, storage_scale=run.storage_scale,
                                            relaxation_delay_s=run.relaxation_delay_s)
            inputs[f"{role}_regions"] = tuple((signal, region) for signal in signals)
            inputs[f"{role}_acquisition"] = acquisition_checks(inputs[f"{role}_raw"], signals)
        assay = link_external_checks(self.assay, **inputs)
        report = self.report(assay)
        self.assertEqual(report["status"], "cannot be assessed")
        self.assertIsNone(report["purity"]["value"])
        self.assertEqual(report["budget"], [])
        self.assertEqual((report["c_s"], report["concentration_unit"]), (original.c_s, "mmol/L"))
        for name in ("c_a", "variance", "u_c", "u_c/w", "U"):
            self.assertIsNone(report[name])
        for name in ("J_a", "J_s", "R"):
            self.assertEqual(report[name], getattr(original, name))
        self.assertEqual(report["prerequisites"], [asdict(check) for check in assay.prerequisites])
        reason = "H146 calibrant region 1: baseline: distorted"
        self.assertEqual(report["reasons"], [reason])
        rendered = render_text(report)
        for line in ("cannot be assessed", reason, "J_a = 0.017390308", "J_s = 0.010959454", "R = 1.586786"):
            self.assertIn(line, rendered)
        self.assertNotIn("Purity =", rendered)
        self.assertIn("c_s = 2.7133857 mmol/L (prepared)", rendered)

    def test_nonfinite_inputs_become_null(self):
        sample = replace(self.assay.sample, area=inf, receiver_gain=-inf, temperature=nan)
        correlation = replace(self.assay.correlations[0], r=-inf)
        assay = replace(self.assay, sample=sample, R_relative_u=nan, N_a=10**400,
                        K=replace(self.assay.K, value=inf, standard_u=nan),
                        correlations=(correlation, None), k=replace(self.assay.k, value=nan))
        report = self.report(assay)
        self.assertEqual(report["status"], "cannot be assessed")
        for name in ("area", "receiver_gain", "temperature"):
            self.assertIsNone(report["sample"][name])
        for name in ("R", "R_relative_u", "N_a", "K"):
            self.assertIsNone(report["inputs"][name]["value"])
            self.assertIsNone(report["inputs"][name]["u"])
        self.assertIsNone(report["covariances"][0]["r"])
        self.assertIsNone(report["k"]["value"])
        self.assertIn("correlations must contain Correlation records", render_text(report))
