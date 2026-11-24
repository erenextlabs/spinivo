import json
import unittest
from dataclasses import asdict, replace
from math import nan
from pathlib import Path

from spinivo.internal_standard import CALCULATION_ONLY, Prerequisite, calculate_internal_standard, link_spectrum_checks
from spinivo.internal_standard_report import purity_report, render_json, render_text
from tests.test_internal_standard import acquisition_checks, clean_regions, load_assay_case


class InternalStandardReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assay, cls.expected = load_assay_case()
        cls.run_id = "FID06"
        raw, clean, _ = clean_regions()
        clean = replace(clean, run=cls.run_id)
        cls.analyte = tuple((name, replace(clean, nominal_area=area))
                            for name, area in zip(cls.assay.analyte_signals, cls.assay.analyte_areas))
        cls.standard = (("G5", replace(clean, nominal_area=cls.assay.standard_areas[0])),)
        raw = replace(raw, source=Path(cls.run_id), pulse_program="single_pulse")
        cls.checks = acquisition_checks(raw, (*cls.assay.analyte_signals, *cls.assay.standard_signals))
        cls.snr = cls.checks[1]

    def link(self, *, analyte=None, checks=None):
        starting = replace(self.assay, analyte_areas=None, standard_areas=None,
                           prerequisites=(Prerequisite("same tube", "pass", "Documented"),))
        return link_spectrum_checks(starting, run=self.run_id, analyte_regions=iter(analyte or self.analyte),
                                    standard_regions=self.standard, acquisition=self.checks if checks is None else checks)

    def test_linked_json_and_text(self):
        unrelated = replace(self.snr, signal="other signal", run="other run", status="flag")
        assay = self.link(checks=(*self.checks, unrelated))
        self.assertEqual((assay.analyte_areas, assay.standard_areas), (self.assay.analyte_areas, self.assay.standard_areas))
        result = calculate_internal_standard(assay)
        self.assertEqual(result, calculate_internal_standard(self.assay))
        report = purity_report(result, assay)
        decoded = json.loads(render_json(report))
        self.assertEqual(decoded, report)
        self.assertEqual(decoded["schema"], "spinivo.internal_standard_purity/2")
        self.assertEqual(decoded["method"], "internal standard, same tube")
        self.assertEqual((decoded["status"], decoded["reasons"]), ("cannot be assessed", [CALCULATION_ONLY]))
        for name in ("assay_id", "analyte_tube_id", "standard_tube_id"):
            self.assertEqual(decoded[name], getattr(assay, name))
        self.assertEqual(decoded["purity"]["unit"], "fraction")
        self.assertAlmostEqual(decoded["purity"]["value"], self.expected["purity"], places=10)
        for field, expected in (("I_a", "I_a"), ("I_s", "I_s"), ("R", "R"), ("u_c", "standard_u"),
                                ("u_c/w", "relative_standard_u"), ("U", "expanded_u"), ("variance", "variance")):
            self.assertAlmostEqual(decoded[field], self.expected[expected], places=8)
        for row in decoded["budget"]:
            name = row["name"]
            self.assertAlmostEqual(row["sensitivity"], self.expected["sensitivities"][name], places=8)
            self.assertAlmostEqual(row["contribution"], self.expected["contributions"][name], places=8)
        for name, sources in (("analyte_areas", ["FID06 H5 analyte region 1", "FID06 I5 analyte region 2"]),
                              ("standard_areas", ["FID06 G5 standard region 1"])):
            self.assertEqual(decoded["inputs"][name], {"value": list(getattr(self.assay, name)), "u": None,
                                                       "unit": "intensity*Hz", "source": sources})
        self.assertEqual(decoded["inputs"]["N_a"]["value"], assay.N_a)
        self.assertEqual(decoded["inputs"]["m_s"]["source_ids"], ["balance:fixture-FID06"])
        self.assertEqual(decoded["inputs"]["w_s"], {"value": 0.999, "u": 0.0005, "unit": "fraction",
                                                    "source": assay.w_s.source})
        covariance, = decoded["covariances"]
        self.assertEqual((covariance["inputs"], covariance["r"], covariance["source"]),
                         (["m_s", "m_a"], 0.30, assay.correlations[0].source))
        self.assertAlmostEqual(covariance["term"], self.expected["covariance_term"], places=14)
        self.assertEqual(decoded["k"], asdict(assay.k))
        self.assertEqual(decoded["prerequisites"], [asdict(check) for check in assay.prerequisites])
        self.assertEqual(len(assay.prerequisites), 13)
        self.assertTrue(all(check.status == "pass" for check in assay.prerequisites))
        rendered = render_text(decoded)
        for line in ("cannot be assessed", CALCULATION_ONLY, *self.expected["text"], "Contribution", "w_s",
                     "Variance term", "m_s, m_a"):
            self.assertIn(line, rendered)

    def test_linked_blockers_keep_areas(self):
        signal, region = self.analyte[0]
        flagged = replace(region, status="flag", flags=(("baseline", "distorted"),), gaps=(("capture", "no evidence"),))
        wrong_run = replace(region, run="other run")
        unknown_run = replace(region, run=None)
        for assay, status, reason in (
            (self.link(analyte=((signal, flagged), self.analyte[1])), "flag", "baseline: distorted. capture: no evidence"),
            (self.link(analyte=((signal, wrong_run), self.analyte[1])), "cannot be assessed",
             "H5 analyte region 1: Run mismatch: expected 'FID06', got 'other run'"),
            (self.link(analyte=((signal, unknown_run), self.analyte[1])), "cannot be assessed",
             "H5 analyte region 1: region run unknown"),
            (self.link(checks=self.checks[1:]), "cannot be assessed", "H5 recovery: Missing recovery check"),
            (self.link(checks=(*self.checks, replace(self.snr, run="wrong"))), "cannot be assessed",
             "Run mismatch: expected 'FID06', got 'wrong'"),
            (self.link(checks=(*self.checks, replace(self.snr, status="flag", reason="Low S/N"))), "flag", "Low S/N"),
        ):
            with self.subTest(reason=reason):
                result = calculate_internal_standard(assay)
                report = json.loads(render_json(purity_report(result, assay)))
                self.assertEqual(report["status"], "cannot be assessed")
                self.assertIsNone(report["purity"]["value"])
                self.assertIsNone(report["u_c"])
                self.assertTrue(any(check.status == status for check in assay.prerequisites))
                for field in ("I_a", "I_s", "R"):
                    self.assertAlmostEqual(report[field], self.expected[field], places=8)
                rendered = render_text(report)
                for value in ("cannot be assessed", reason, "I_a = 3.6047681", "I_s = 1.2346185", "R = 2.9197424"):
                    self.assertIn(value, rendered)
                self.assertNotIn(" %", rendered)

    def test_mixed_linked_area_units_block_purity(self):
        (signal, region), = self.standard
        different_unit = replace(region, area=replace(region.area, unit="intensity*ppm"))
        assay = link_spectrum_checks(self.assay, run=self.run_id, analyte_regions=self.analyte,
                                     standard_regions=((signal, different_unit),), acquisition=self.checks)
        result = calculate_internal_standard(assay)
        reason = ("Region units differ (H5 analyte region 1 = 'intensity*Hz', I5 analyte region 2 = 'intensity*Hz', "
                  "G5 standard region 1 = 'intensity*ppm')")
        self.assertIn(Prerequisite("linked area units", "cannot be assessed", reason), assay.prerequisites)
        self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
        self.assertIn(f"linked area units: {reason}", result.reasons)
        self.assertIsNone(assay.area_unit)

    def test_mismatched_area_sources_align_in_report(self):
        assay = replace(self.assay, analyte_area_sources=("H5", "I5", "extra"), standard_area_sources=())
        result = calculate_internal_standard(assay)
        report = json.loads(render_json(purity_report(result, assay)))
        self.assertEqual(report["status"], "cannot be assessed")
        self.assertIsNone(report["purity"]["value"])
        self.assertEqual(report["inputs"]["analyte_areas"]["source"], ["H5", "I5"])
        self.assertEqual(report["inputs"]["standard_areas"]["source"], [None])
        self.assertEqual(report["reasons"], list(result.reasons))
        for reason in report["reasons"]:
            self.assertIn(reason, render_text(report))

    def test_negative_ratio_has_nonnegative_uncertainty(self):
        assay = replace(self.assay, analyte_areas=(-3.0,))
        result = calculate_internal_standard(assay)
        report = purity_report(result, assay)
        self.assertEqual(report["status"], "cannot be assessed")
        self.assertLess(report["R"], 0)
        self.assertEqual(report["inputs"]["R"]["u"], abs(result.R) * assay.R_relative_u)

    def test_missing_coverage_and_nonfinite_inputs(self):
        assay = replace(self.assay, k=None)
        report = purity_report(calculate_internal_standard(assay), assay)
        self.assertIsNone(report["U"])
        self.assertIsNone(report["inputs"]["analyte_areas"]["unit"])
        self.assertIsNone(report["inputs"]["analyte_areas"]["source"])
        self.assertIn("coverage factor k was not supplied", render_text(report))
        assay = replace(self.assay, analyte_areas=(nan,), N_a=10**400)
        report = json.loads(render_json(purity_report(calculate_internal_standard(assay), assay)))
        self.assertEqual(report["inputs"]["analyte_areas"]["value"], [None])
        self.assertIsNone(report["inputs"]["N_a"]["value"])
        for reason in report["reasons"]:
            self.assertIn(reason, render_text(report))

    def test_nonstring_units_and_sources_keep_blocking_reasons(self):
        correlation, = self.assay.correlations
        assay = replace(self.assay, w_s=replace(self.assay.w_s, unit=b"%", source=1),
                        k=replace(self.assay.k, source=b"coverage"), R_u_source=b"repeatability",
                        correlations=(replace(correlation, source=b"correlation"),),
                        area_unit=b"intensity*Hz", analyte_area_sources=(b"H5", "I5"))
        result = calculate_internal_standard(assay)
        report = json.loads(render_json(purity_report(result, assay)))
        self.assertEqual(report["status"], "cannot be assessed")
        self.assertIsNone(report["inputs"]["w_s"]["unit"])
        self.assertIsNone(report["inputs"]["w_s"]["source"])
        self.assertIsNone(report["inputs"]["R"]["source"])
        self.assertIsNone(report["inputs"]["analyte_areas"]["unit"])
        self.assertEqual(report["inputs"]["analyte_areas"]["source"], [None, "I5"])
        self.assertIsNone(report["covariances"][0]["source"])
        self.assertIsNone(report["k"]["source"])
        self.assertTrue(any("w_s.unit" in reason for reason in report["reasons"]))
        self.assertEqual(report["reasons"], list(result.reasons))
        for reason in result.reasons:
            self.assertIn(reason, render_text(report))
