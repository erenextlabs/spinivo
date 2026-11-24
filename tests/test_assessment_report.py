import json
import unittest
from dataclasses import asdict, fields, replace
from math import nan
from pathlib import Path
from unittest.mock import patch

from spinivo import __version__, assess_record, load_record
from spinivo.acquisition import PILOT_THRESHOLDS
from spinivo.assessment_report import NOT_ASSESSED, assessment_report, render_json, render_text
from spinivo.external_calibration import TRANSFER_TOLERANCES
from spinivo.integral_checks import PILOT
from spinivo.integral_components import check_capture_pair
from spinivo.internal_standard import calculate_internal_standard
from spinivo.pipeline import RunEvidence
from spinivo.processing import QUANTITATIVE_SETTINGS, process
from tests import test_pipeline as spectral
from tests.test_assay_record import ENTERED, RECORDS


class AssessmentReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.internal, cls.external = spectral.PipelineTests, spectral.ExternalPipelineTests
        cls.internal.setUpClass()
        cls.external.setUpClass()
        cls.record = load_record(RECORDS / f"{ENTERED}.json")
        cls.entered = assess_record(cls.record)

    def report(self, obj, record=None):
        report = assessment_report(obj, record)
        decoded = json.loads(render_json(report), parse_constant=self.fail)
        self.assertEqual(decoded, report)
        return decoded, render_text(decoded)

    def test_strict_json_and_text_agree_for_every_mode(self):
        for obj, record, method, mode, status in (
            (self.internal.complete, None, "internal_standard", "internal_standard", "pass"),
            (self.external.complete, None, "external_calibration", "external_calibration", "pass"),
            (self.entered, self.record, "internal_standard", "calculation_only", "cannot be assessed"),
        ):
            with self.subTest(mode=mode):
                report, text = self.report(obj, record)
                method_report = report["method_report"]
                self.assertEqual((report["schema"], report["spinivo_version"], report["method"], report["mode"]),
                                 ("spinivo.assessment/1", __version__, method, mode))
                self.assertEqual((report["assay_id"], report["status"], report["reasons"]),
                                 (method_report["assay_id"], status, method_report["reasons"]))
                self.assertIn(f'{report["assay_id"]} ({method}, {mode} mode): {status}', text)
                self.assertIn(f'R = {method_report["R"]:.8g}', text)
                for reason in report["reasons"]:
                    self.assertIn(f"    {reason}", text)
                for check in method_report["prerequisites"]:
                    self.assertIn(f'{check["name"]}: {check["status"]}: {check["reason"]}', text)
                if mode == "calculation_only":
                    self.assertEqual((report["profiles"], report["runs"], report["not_assessed"], report["purity"]),
                                     (None, None, list(NOT_ASSESSED), {"value": None, "unit": "fraction"}))
                    self.assertIn("Not assessed in calculation_only mode\n    " + "\n    ".join(NOT_ASSESSED), text)
                    self.assertNotIn("baseline_shape", report)
                    self.assertNotIn("Baseline shape (", text)
                    continue
                shape = report["baseline_shape"]
                self.assertEqual((shape["method_version"], shape["status"], shape["reason"]),
                                 ("baseline-shape/1", "computed", None))
                value = shape["relative_standard_uncertainty_percent"]
                self.assertIn(f"Baseline shape (baseline-shape/1, diagnostic): computed, {value:.8g} %", text)
                self.assertEqual(text.count("Baseline shape ("), 1)
                for entry in shape["regions"]:
                    (region,) = (item for item in report["runs"][entry["run"]]["regions"]
                                 if item["role"] == entry["role"] and item["signal"] == entry["signal"])
                    self.assertEqual(set(entry["candidate_areas"]), {"production", "anchor_median", "bracketing_line"})
                    self.assertEqual(entry["candidate_areas"]["production"], region["report"]["nominal_area"]["value"])
                self.assertEqual(report["purity"], {"value": method_report["purity"]["value"], "unit": "fraction"})
                self.assertIn(f'Checked purity: {report["purity"]["value"]:.8g} fraction', text)
                self.assertEqual(report["profiles"]["baseline"],
                                 {"algorithm_version": "locked-anchors/1", "limit": 0.017})
                for run in report["runs"].values():
                    self.assertEqual(run["processing"]["autophase_method"], "lineshape")
                    for check in run["acquisition"]:
                        self.assertIn(f'{check["signal"]} {check["check"]}: {check["status"]}, value = '
                                      f'{check["value"]:.8g} {check["unit"]}, limit = {check["limit"]:.8g}: '
                                      f'{check["reason"]}', text)
                        self.assertTrue(all(f'{item["name"]} = ' in text and f'({item["source"]})' in text
                                            for item in check["inputs"]))
                    for region in run["regions"]:
                        area = region["report"]["nominal_area"]
                        self.assertIn(f'Area: {area["value"]:.8g} {area["unit"]}', text)

    def test_baseline_shape_leaves_the_existing_report_bit_identical(self):
        for assay in (self.internal, self.external):
            with self.subTest(assay=assay.__name__):
                computed = assessment_report(assay.complete)
                with patch("spinivo.processing._baseline_shape_areas", return_value=None):
                    unavailable = assessment_report(assay.assess())
                self.assertEqual(computed.pop("baseline_shape")["status"], "computed")
                self.assertEqual(unavailable.pop("baseline_shape")["status"], "cannot be assessed")
                self.assertEqual(render_json(computed), render_json(unavailable))

    def test_processing_reports_fitted_phase_and_pivot(self):
        raw = replace(self.internal.raw, fid=self.internal.raw.fid * 1j)
        report, _ = self.report(self.internal.assess(raw))
        bounds = tuple(signal.bounds_ppm for signal in self.internal.signals.values())
        applied = process(raw, QUANTITATIVE_SETTINGS, regions_ppm=bounds).settings
        actual = report["runs"]["run"]["processing"]
        self.assertEqual(tuple(actual[key] for key in ("pivot_ppm", "p0_deg", "p1_deg")),
                         (applied.pivot_ppm, applied.p0_deg, applied.p1_deg))
        self.assertNotEqual(applied.p0_deg, QUANTITATIVE_SETTINGS.p0_deg)

    def test_integral_version_and_algorithm_stay_as_evaluated(self):
        assessment = self.internal.complete
        with patch("spinivo.integral_checks.PILOT", replace(PILOT, version="later")), \
                patch("spinivo.pipeline.QUANTITATIVE_SETTINGS", replace(QUANTITATIVE_SETTINGS, lock_spike_anchors=False)):
            report, _ = self.report(assessment)
        self.assertEqual(report["profiles"]["baseline"]["algorithm_version"], "locked-anchors/1")
        profile = report["profiles"]["integral"]
        self.assertEqual(profile["version"], "pilot-3")
        for entry in report["runs"]["run"]["regions"]:
            self.assertEqual(entry["report"]["thresholds"], profile)
            self.assertEqual({obs["threshold_version"] for obs in entry["report"]["observations"]}, {"pilot-3"})

    def test_baseline_recovery_uses_applied_limit(self):
        limit = 1e-12
        with patch("spinivo.pipeline._BASELINE_SENSITIVITY_LIMIT", limit):
            assessment = self.internal.assess()
        report, _ = self.report(assessment)
        self.assertEqual(report["profiles"]["baseline"]["limit"], limit)
        recovery = assessment.runs[0].regions[0][2].baseline.details.recovery
        self.assertEqual((recovery.limit, recovery.flag), (limit, "baseline_sensitive"))

    def test_external_tolerances_are_captured_at_evaluation(self):
        external = self.external
        sample = replace(external.raws["sample"], temperature_k=298.35)
        assay = replace(external.assay, sample=replace(external.assay.sample, temperature=25.2))
        with patch.dict(TRANSFER_TOLERANCES, temperature=0.25):
            assessment = external.assess(sample_fid=sample, assay=assay)
        with patch.dict(TRANSFER_TOLERANCES, temperature=0.1):
            report, _ = self.report(assessment)
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["method_report"]["tolerances"], {"temperature": 0.25, "observe_frequency_mhz": 1e-4})

    def test_calculation_only_requires_matching_result_and_labels_value(self):
        assay = replace(self.record.assay, analyte_areas=(2.0, 3.0))
        with self.assertRaisesRegex(ValueError, "does not match"):
            self.report(self.entered, replace(self.record, assay=assay))
        report, text = self.report(self.entered, self.record)
        self.assertIsNone(report["purity"]["value"])
        self.assertEqual(report["method_report"]["purity"]["value"], self.entered.value)
        self.assertRegex(text, r"Purity = [0-9.]+ %, u_c = [0-9.]+ % \(calculation only, not checked\)")

    def test_text_keeps_all_observation_values_and_reasons(self):
        signal = self.internal.signals["H5"]
        evidence = replace(signal.evidence, target_source=None)
        for result in (self.internal.complete, self.internal.assess(
                analyte=(replace(signal, evidence=evidence), self.internal.signals["I5"]))):
            report, text = self.report(result)
            for entry in report["runs"]["run"]["regions"]:
                for obs in entry["report"]["observations"]:
                    value, limit = ("missing" if obs[key] is None else f"{obs[key]:.8g}" for key in ("value", "limit"))
                    self.assertIn(f'{obs["check"]}: {obs["status"]}, value = {value} {obs["unit"]}, '
                                  f'limit = {limit}: {obs["reason"]}', text)

    def test_capture_pairs_keep_both_signals_and_observation(self):
        for assessment in (self.internal.complete, self.external.complete):
            report, text = self.report(assessment)
            entries = [entry for run in assessment.runs for entry in run.regions]
            first, second = entries[0], entries[-1]
            expected = check_capture_pair(first[2], second[2], assessment.profiles.integral)
            pair = report["capture_pairs"][0]
            self.assertEqual(pair["signals"], [first[1], second[1]])
            for key in ("check", "status", "value", "unit", "limit", "threshold_version", "reason"):
                self.assertEqual(pair[key], getattr(expected, key))
            self.assertIn(f'{first[1]}/{second[1]} capture_pair: {expected.status}, value = {expected.value:.8g} '
                          f'{expected.unit}, limit = {expected.limit:.8g}: {expected.reason}', text)

    def test_report_keeps_the_profile_applied_at_evaluation(self):
        applied = replace(PILOT, unexpected_peak_sd=999.0, min_noise_points=64)
        with patch("spinivo.pipeline.PILOT", applied):
            assessment = self.internal.assess()
        report, text = self.report(assessment)
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["profiles"]["integral"], asdict(replace(applied, version="pilot-3")))
        for region in report["runs"]["run"]["regions"]:
            profile = region["report"]["thresholds"]
            self.assertEqual((profile["unexpected_peak_sd"], profile["min_noise_points"]), (999, 64))
        self.assertIn("min_noise_points = 64", text)
        self.assertIn("unexpected_peak_sd = 999", text)
        for override in (applied, assessment.profiles.integral):
            self.assertEqual(assessment_report(assessment, thresholds=override), report)
        for override in (PILOT, replace(applied, min_noise_points=32),
                         replace(applied, max_pair_capture_mismatch=0.01)):
            with self.subTest(override=override), self.assertRaisesRegex(ValueError, "does not match"):
                assessment_report(assessment, thresholds=override)
        with self.assertRaisesRegex(ValueError, "does not match"):
            assessment_report(self.entered, self.record, thresholds=PILOT)
        default = self.report(self.internal.complete)[0]["profiles"]["integral"]
        self.assertEqual((default["unexpected_peak_sd"], default["min_noise_points"]), (5, 32))

    def test_missing_run_id_keeps_the_independent_resolution_flag(self):
        external = self.external
        assay = replace(external.assay, sample=replace(external.assay.sample, run_id=None))
        strict = replace(PILOT_THRESHOLDS, min_points_per_fwhm=1000.0)
        with patch("spinivo.pipeline.PILOT_THRESHOLDS", strict):
            assessment = external.assess(assay=assay)
        report, text = self.report(assessment)
        self.assertEqual((report["status"], report["purity"]["value"]), ("cannot be assessed", None))
        self.assertIn("sample run_id is missing", report["reasons"])
        self.assertEqual(report["profiles"]["acquisition"], asdict(strict))
        name, = external.assay.sample_signals
        linked = {check["name"] for check in report["method_report"]["prerequisites"]}
        self.assertNotIn(f"sample {name} resolution", linked)
        resolution, = (check for check in report["runs"]["sample"]["acquisition"] if check["check"] == "resolution")
        self.assertEqual((resolution["signal"], resolution["status"], resolution["limit"]), (name, "flag", 1000))
        self.assertIn(f"{name} resolution: flag, value = {resolution['value']:.8g} points / FWHM, limit = 1000", text)

    def test_direct_calls_take_run_ids_from_the_assay_unless_given(self):
        internal, external = self.internal.assay, self.external.assay
        self.assertEqual(internal.analyte_tube_id, internal.standard_tube_id)
        self.assertNotIn(None, (internal.analyte_tube_id, *(getattr(getattr(external, role), name)
                                for role in ("sample", "calibrant") for name in ("run_id", "preparation_id"))))
        for assessment, expected in (
            (self.internal.complete, {"run": ("FID06", internal.analyte_tube_id)}),
            (self.external.complete, {role: (getattr(external, role).run_id, getattr(external, role).preparation_id)
                                      for role in ("sample", "calibrant")}),
            (self.internal.assess(run_id="tube-run", preparation_id=None), {"run": ("tube-run", None)}),
            (self.external.assess(sample_run_id=None, calibrant_preparation_id="other-prep"),
             {"sample": (None, external.sample.preparation_id), "calibrant": (external.calibrant.run_id, "other-prep")}),
        ):
            report, text = self.report(assessment)
            for label, ids in expected.items():
                run = report["runs"][label]
                self.assertEqual((run["run_id"], run["preparation_id"]), ids)
                self.assertIn(f'Run {label}: {run["source"]}, run_id = {ids[0] or "missing"}, '
                              f'preparation_id = {ids[1] or "missing"}', text)

    def test_run_ids_follow_the_existing_run_fields(self):
        names = [field.name for field in fields(RunEvidence)]
        self.assertEqual(names, ["run", "source", "settings", "regions", "acquisition", "run_id", "preparation_id"])
        self.assertEqual(RunEvidence("run", "file", QUANTITATIVE_SETTINGS).settings, QUANTITATIVE_SETTINGS)

    def test_text_renders_runs_without_id_keys(self):
        report, _ = self.report(self.internal.complete)
        for run in report["runs"].values():
            del run["run_id"], run["preparation_id"]
        self.assertIn(f'Run run: {report["runs"]["run"]["source"]}, run_id = missing, preparation_id = missing',
                      render_text(report))

    def test_withheld_purity_and_missing_evidence_stay_visible(self):
        blocked = self.internal.assess(recovery=replace(spectral.RECOVERY, calibration_source=None))
        report, text = self.report(blocked)
        self.assertEqual((report["status"], report["purity"]["value"], report["method_report"]["purity"]["value"]),
                         ("cannot be assessed", None, None))
        self.assertIn("Checked purity: cannot be assessed", text)
        self.assertNotIn("Purity =", text)
        self.assertEqual(report["method_report"]["R"], self.internal.complete.result.R)
        with patch("spinivo.pipeline.read", side_effect=OSError("corrupt data")):
            unread, _ = self.report(self.internal.assess(Path("unreadable.jdf"), run_id=None))
        tube = self.internal.assay.analyte_tube_id
        self.assertEqual(unread["runs"], {"run": {"source": None, "run_id": None, "preparation_id": tube,
                                                  "processing": None, "acquisition": None, "regions": None}})
        h5, i5 = (self.internal.signals[name] for name in ("H5", "I5"))
        unbounded, text = self.report(self.internal.assess(analyte=(replace(h5, bounds_ppm=None), i5)))
        self.assertEqual(unbounded["runs"]["run"]["regions"][0], {"role": "analyte", "signal": "H5", "report": None})
        self.assertIn("H5 analyte: no region result", text)
        assay = replace(self.record.assay, analyte_areas=(nan, 1.0))
        entered, _ = self.report(calculate_internal_standard(assay), replace(self.record, assay=assay))
        self.assertEqual(entered["method_report"]["inputs"]["analyte_areas"]["value"], [None, 1.0])
        with self.assertRaisesRegex(ValueError, "does not match"):
            assessment_report(self.internal.complete, self.record)
        with self.assertRaisesRegex(TypeError, "requires its record"):
            assessment_report(self.internal.complete.result)
