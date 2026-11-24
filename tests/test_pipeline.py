import unittest
from contextlib import chdir
from dataclasses import asdict, dataclass, fields, replace
from math import hypot, log, sqrt
from os import environ
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import DEFAULT, patch

import numpy as np

from spinivo import DeclaredSignal, UnsupportedInput, assess_external_calibration, assess_internal_standard
from spinivo.acquisition import RecoveryEvidence
from spinivo.external_calibration import ExternalRun, calculate_external_calibration
from spinivo.external_calibration_report import purity_report as external_report, render_text as external_text
from spinivo.integral_components import ComponentEvidence, Line, check_region
from spinivo.internal_standard import CALCULATION_ONLY, Prerequisite, calculate_internal_standard, link_spectrum_checks
from spinivo.internal_standard_report import purity_report, render_text
from spinivo.measurement import measure_region
from spinivo.processing import QUANTITATIVE_SETTINGS, BaselineShapeAreas, process
from tests.integral_fixtures import ANCHORS, CARRIER_PPM, NOISE_BOUNDS_HZ, OBSERVE_MHZ, _fid, make_raw, ppm_bounds
from tests.test_external_calibration import load_assay_case as load_external_case
from tests.test_internal_standard import acquisition_checks, clean_regions, load_assay_case
from tests.test_readers import FIXTURES, JEOL, JEOL_ARRAY
from tools.compare_references import load_references

# FID06 signals as ideal 1 Hz lines, by center in Hz from the carrier, each integrated over +/-50 Hz.
CENTERS = {"H5": -1200, "I5": 0, "G5": 1200}
RECOVERY = RecoveryEvidence(calibrated_flip_angle_deg=90.0, calibration_source="pulse log")


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assay, cls.expected = load_assay_case()
        areas = dict(zip(CENTERS, (*cls.assay.analyte_areas, *cls.assay.standard_areas)))
        # Noise-free lines: the noise window holds only their wings, so S/N is far above the limit.
        fid = _fid((CENTERS[name], area, 1) for name, area in areas.items())
        cls.raw = replace(make_raw({"id": "FID06"}, ideal=True), fid=fid, relaxation_delay_s=60.0,
                          pulse_width_us=10.0, pulse_program="single_pulse")
        cls.signals = {}
        for name, area in areas.items():
            center = CENTERS[name]
            evidence = ComponentEvidence((Line(CARRIER_PPM + center / OBSERVE_MHZ, area, 1.0),), "fixture",
                                         satellite_convention="absent", satellite_source="fixture")
            cls.signals[name] = DeclaredSignal(name, ppm_bounds((center - 50, center + 50)), evidence, 2.0, "T1 log")
        cls.complete = cls.assess()

    @classmethod
    def assess(cls, fid=None, assay=None, analyte=None, standard=None, **changes):
        return assess_internal_standard(
            cls.raw if fid is None else fid, cls.assay if assay is None else assay,
            analyte=(cls.signals["H5"], cls.signals["I5"]) if analyte is None else analyte,
            standard=(cls.signals["G5"],) if standard is None else standard,
            **{"noise_bounds_ppm": ppm_bounds(NOISE_BOUNDS_HZ), "recovery": RECOVERY, "run_id": "FID06", **changes})

    def assertWithheld(self, assessment, *reasons):
        result = assessment.result
        self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
        for reason in reasons:
            self.assertIn(reason, result.reasons)
        return result

    def test_complete_evidence_passes_with_hand_worked_purity(self):
        assay, result = self.complete.assay, self.complete.result
        self.assertEqual((result.status, result.reasons), ("pass", ()))
        self.assertEqual(len(assay.prerequisites), 14)
        self.assertTrue(all(check.status == "pass" for check in assay.prerequisites))
        self.assertAlmostEqual(result.R, self.expected["R"], delta=1e-4 * self.expected["R"])
        self.assertAlmostEqual(result.value, self.expected["purity"], delta=1e-4)
        self.assertNotEqual(assay.analyte_areas, self.assay.analyte_areas)
        self.assertEqual(assay.analyte_area_sources, ("FID06 H5 analyte region 1", "FID06 I5 analyte region 2"))
        report = purity_report(self.complete)
        self.assertEqual((report["schema"], report["status"]), ("spinivo.internal_standard_purity/2", "pass"))
        for line in self.expected["text"]:
            self.assertIn(line, render_text(report))
        direct = calculate_internal_standard(assay)
        self.assertEqual((direct.status, direct.reasons), ("cannot be assessed", (CALCULATION_ONLY,)))
        self.assertEqual((direct.value, direct.budget), (result.value, result.budget))

    def test_repeated_assessments_are_identical(self):
        self.assertEqual(self.assess(), self.complete)
        self.assertEqual(self.assess(assay=self.complete.assay), self.complete)

    def test_baseline_shape_sums_signals_before_the_same_tube_ratio(self):
        for standard, departure in ((BaselineShapeAreas(20, 22, 20), 0), (BaselineShapeAreas(20, 20, 20), 0.1)):
            with self.subTest(standard=standard):
                areas = (BaselineShapeAreas(10, 14, 9), BaselineShapeAreas(30, 30, 31), standard)
                with patch("spinivo.pipeline.baseline_shape_areas", side_effect=areas):
                    assessed = self.assess()
                shape = assessed.baseline_shape
                self.assertEqual(shape.status, "computed")
                self.assertAlmostEqual(shape.relative_standard_uncertainty_percent, 100 * departure / sqrt(3),
                                       places=12)
                self.assertEqual(purity_report(assessed), purity_report(self.complete))
                self.assertEqual(assessed.result, self.complete.result)

    def test_invalid_baseline_shape_candidates_are_diagnostic_only(self):
        for area, reason in ((None, "missing accepted baseline anchor"), (0.0, "finite and positive"),
                             (float("inf"), "finite and positive"), (float("nan"), "finite and positive")):
            with self.subTest(area=area):
                with patch("spinivo.pipeline.baseline_shape_areas", return_value=BaselineShapeAreas(10, 10, area)):
                    assessed = self.assess()
                shape = assessed.baseline_shape
                self.assertEqual(shape.status, "cannot be assessed")
                self.assertIn(reason, shape.reason)
                self.assertIsNone(shape.relative_standard_uncertainty_percent)
                self.assertEqual(assessed.result, self.complete.result)

    def test_baseline_shape_requires_quantitative_settings(self):
        with patch("spinivo.pipeline.QUANTITATIVE_SETTINGS", replace(QUANTITATIVE_SETTINGS, lock_spike_anchors=False)):
            shape = self.assess().baseline_shape
        self.assertEqual(shape.status, "not computed")
        self.assertIn("QUANTITATIVE_SETTINGS", shape.reason)
        self.assertEqual(shape.regions, ())

    def test_baseline_shape_requires_a_source_for_a_nonzero_reference_shift(self):
        for source, status in ((None, "not computed"), ("  ", "not computed"), ("fixture", "computed")):
            with self.subTest(source=source):
                settings = replace(QUANTITATIVE_SETTINGS, reference_shift_ppm=0.01, reference_source=source)
                with patch("spinivo.pipeline.QUANTITATIVE_SETTINGS", settings):
                    shape = self.assess().baseline_shape
                self.assertEqual(shape.status, status)

    def test_internal_pipeline_requires_a_run_id(self):
        inputs = {"analyte": (self.signals["H5"], self.signals["I5"]), "standard": (self.signals["G5"],),
                  "noise_bounds_ppm": ppm_bounds(NOISE_BOUNDS_HZ), "recovery": RECOVERY}
        for ids in ({}, {"run_id": None}, {"run_id": "  "}, {"run_id": "NOT FOUND"}):
            with self.subTest(ids=ids):
                blocked = assess_internal_standard(self.raw, self.assay, **inputs, **ids)
                result = self.assertWithheld(blocked, "run_id is missing")
                self.assertEqual(result.reasons, ("run_id is missing",))
                self.assertEqual(calculate_internal_standard(blocked.assay).reasons, result.reasons)
                self.assertEqual((result.I_a, result.I_s, result.R),
                                 (self.complete.result.I_a, self.complete.result.I_s, self.complete.result.R))
                self.assertEqual(blocked.runs[0].regions, self.complete.runs[0].regions)

    def test_signals_are_measured_in_assay_declaration_order(self):
        for analyte, standard in ((("H5", "I5"), ("G5",)), (("H5",), ("I5", "G5"))):
            with self.subTest(analyte=analyte, standard=standard):
                assay = replace(self.assay, analyte_signals=analyte, standard_signals=standard)
                a, s = (tuple(self.signals[name] for name in names) for names in (analyte, standard))
                expected = self.assess(assay=assay, analyte=a, standard=s)
                reversed_inputs = self.assess(assay=assay, analyte=reversed(a), standard=reversed(s))
                self.assertEqual(expected.result.status, "pass")
                self.assertEqual(reversed_inputs, expected)

    def test_supplied_signal_names_must_match_the_assay_by_role(self):
        h5, i5, g5 = (self.signals[name] for name in ("H5", "I5", "G5"))
        for analyte, standard, reasons in (
            ((h5,), (g5,), ("I5 analyte signal: Missing declared signal",)),
            ((h5, i5), (), ("G5 standard signal: Missing declared signal",)),
            ((h5, i5, replace(g5, name="J5", bounds_ppm=None)), (g5,), ("J5 analyte signal: Undeclared signal",)),
            ((h5, replace(i5, name="J5")), (g5,),
             ("I5 analyte signal: Missing declared signal", "J5 analyte signal: Undeclared signal")),
            ((h5, g5), (i5,), ("G5 analyte signal: Undeclared signal", "I5 standard signal: Undeclared signal")),
        ):
            with self.subTest(reasons=reasons):
                assessment = self.assess(analyte=analyte, standard=standard)
                self.assertWithheld(assessment, *reasons)
                self.assertEqual(assessment.assay.analyte_signals, self.assay.analyte_signals)
                self.assertEqual(assessment.assay.standard_signals, self.assay.standard_signals)

    def test_checked_report_uses_the_assessment_snapshot(self):
        assessment = self.complete
        self.assertEqual(assessment.result.status, "pass")
        edited = replace(assessment.assay, m_a=replace(self.assay.m_a, value=2 * self.assay.m_a.value))
        with self.assertRaisesRegex(ValueError, "checked pass report"):
            purity_report(assessment.result, edited)
        with self.assertRaises(TypeError):
            purity_report(assessment, edited)
        with self.assertRaises(TypeError):
            replace(assessment, _assay=edited)
        report = purity_report(assessment)
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["inputs"]["m_a"]["value"], self.assay.m_a.value)
        self.assertEqual(len(report["covariances"]), len(self.assay.correlations))

    def test_assessment_freezes_caller_collections(self):
        correlation, = self.assay.correlations
        inputs = list(correlation.inputs)
        correlations = [replace(correlation, inputs=inputs)]
        source_ids = list(self.assay.m_s.source_ids)
        prerequisites = [Prerequisite("same tube", "pass", "Documented")]
        names = list(self.assay.analyte_signals)
        assessment = self.assess(assay=replace(self.assay, correlations=correlations, prerequisites=prerequisites,
                                             analyte_signals=names, m_s=replace(self.assay.m_s, source_ids=source_ids)))
        before = purity_report(assessment)
        self.assertEqual(before["status"], "pass")

        inputs.clear()
        correlations.clear()
        source_ids.clear()
        prerequisites.clear()
        names.reverse()
        self.assertEqual(purity_report(assessment), before)

        snapshot = assessment.assay
        for items in (snapshot.correlations, snapshot.correlations[0].inputs, snapshot.m_s.source_ids,
                      snapshot.prerequisites, snapshot.analyte_signals, snapshot.analyte_area_sources,
                      assessment.result.budget.covariances[0].correlation.inputs):
            self.assertIsInstance(items, tuple)

    def test_review_exploits_never_pass(self):
        h5, i5, g5 = (self.signals[name] for name in ("H5", "I5", "G5"))
        for analyte, standard, reason in (
            ((h5, replace(h5, name="I5")), (g5,), "I5 analyte region: Overlaps the H5 analyte region"),
            ((h5, replace(i5, name="H5")), (g5,), "H5 analyte signal: Signal name is declared more than once"),
            ((h5, i5), (replace(g5, name="H5"),), "H5 standard signal: Signal name is declared more than once"),
            ((h5, replace(i5, bounds_ppm=None)), (g5,), "I5 analyte region: Missing region"),
            ((h5, replace(i5, bounds_ppm=(3.0, 3.0))), (g5,), "I5 analyte bounds: Declare two distinct finite ppm bounds"),
        ):
            with self.subTest(reason=reason):
                self.assertWithheld(self.assess(analyte=analyte, standard=standard), reason)

        linked = self.complete.assay
        removed = calculate_internal_standard(replace(linked, analyte_areas=linked.analyte_areas[:1],
                                                      analyte_area_sources=linked.analyte_area_sources[:1]))
        self.assertEqual((removed.status, removed.reasons), ("cannot be assessed", (CALCULATION_ONLY,)))
        self.assertAlmostEqual(removed.value, self.complete.result.value / 2, delta=0.01)
        raw, region, _ = clean_regions()
        region = replace(region, run="FID06")
        raw = replace(raw, source=Path("FID06"), pulse_program="single_pulse")
        doubled = link_spectrum_checks(self.assay, run="FID06", analyte_regions=(("H5", region), ("I5", region)),
                                       standard_regions=(("G5", region),),
                                       acquisition=acquisition_checks(raw, ("H5", "I5", "G5")))
        self.assertTrue(all(check.status == "pass" for check in doubled.prerequisites))
        self.assertEqual(calculate_internal_standard(doubled).reasons, (CALCULATION_ONLY,))

    def test_unsupported_input_propagates(self):
        for fid, reason in ((JEOL_ARRAY, "arrayed_or_2d"), (Path("missing.jdf"), "not_found")):
            with self.subTest(reason=reason):
                with self.assertRaises(UnsupportedInput) as caught:
                    self.assess(fid)
                self.assertEqual(caught.exception.reason, reason)

    def test_reader_and_processing_failures_are_diagnosed(self):
        for error in (OSError, ValueError):
            with self.subTest(error=error), patch("spinivo.pipeline.read", side_effect=error("corrupt data")):
                self.assertWithheld(self.assess(JEOL), "reader: corrupt data")
        outside = replace(self.signals["H5"], bounds_ppm=(50.0, 51.0))
        self.assertWithheld(self.assess(analyte=(outside, self.signals["I5"])),
                            "processing: ppm bounds must span a nonzero interval within the spectral width")

    def test_raw_fid_requires_proton_nucleus(self):
        self.assertWithheld(self.assess(replace(self.raw, nucleus="13C")), "nucleus: Expected 1H, got '13C'")

    def test_placeholder_raw_source_blocks_before_measuring(self):
        with patch("spinivo.pipeline.measure_region") as measure:
            self.assertWithheld(self.assess(replace(self.raw, source="NF")), "source: Missing raw source")
            measure.assert_not_called()

    def test_reference_without_a_source_is_a_named_prerequisite(self):
        self.assertWithheld(self.assess(reference=(0.5,)), "reference: Declare shift_ppm with its source")

    def test_raw_fid_timing_cannot_turn_failed_recovery_into_pass(self):
        short = replace(self.raw, relaxation_delay_s=0.1)
        self.assertWithheld(self.assess(short), "H5 recovery: Outside the pilot limit for the steady-state bound")
        self.assertWithheld(self.assess(replace(short, acquisition_time_s=100.0)),
                            "acquisition_time_s: 100 s does not match complex_points / spectral_width_hz (8 s)")
        for field, value, reason in (
            ("acquisition_time_s", None, "acquisition_time_s: Must be finite and positive"),
            ("spectral_width_hz", 0.0, "spectral_width_hz: Must be finite and positive"),
            ("complex_points", float(self.raw.complex_points), "complex_points: Must be a positive integer"),
        ):
            with self.subTest(field=field):
                self.assertWithheld(self.assess(replace(self.raw, **{field: value})), reason)

    def test_positive_extra_recovery_requires_a_source(self):
        short = replace(self.raw, relaxation_delay_s=0.1)
        extra = replace(RECOVERY, extra_recovery_s=20.0)
        sourced = self.assess(short, recovery=replace(extra, extra_recovery_source="timing log"))
        self.assertEqual(sourced.result.status, "pass")
        independent = tuple(check for check in sourced.assay.prerequisites if not check.name.endswith(" recovery"))
        for source in (None, "  "):
            assessment = self.assess(short, recovery=replace(extra, extra_recovery_source=source))
            result = self.assertWithheld(assessment, "extra_recovery_s: extra_recovery_source is missing")
            self.assertEqual((result.I_a, result.I_s, result.R), (sourced.result.I_a, sourced.result.I_s, sourced.result.R))
            self.assertEqual(assessment.assay.analyte_areas, sourced.assay.analyte_areas)
            self.assertEqual(assessment.assay.standard_areas, sourced.assay.standard_areas)
            checks = assessment.assay.prerequisites
            self.assertIn(Prerequisite("extra_recovery_s", "cannot be assessed", "extra_recovery_source is missing"),
                          checks)
            self.assertTrue(set(independent).issubset(checks))
            self.assertTrue(all(check.status == "cannot be assessed" for check in checks
                                if check.name.endswith(" recovery")))
        self.assertWithheld(self.assess(recovery=replace(extra, extra_recovery_s=True)),
                            "extra_recovery_s: Must be finite and nonnegative")

    def test_recovery_uses_the_sampling_grid_time(self):
        # With a 90 degree pulse the loss is exp(-TR / T1); T1 = 2 s puts the 0.01 limit at TR = 2 ln 100.
        limit = 2 * log(100) - self.raw.acquisition_time_s
        reason = "H5 recovery: Outside the pilot limit for the steady-state bound"
        self.assertEqual(self.assess(replace(self.raw, relaxation_delay_s=limit + 2e-9)).result.status, "pass")
        short = replace(self.raw, relaxation_delay_s=limit - 2e-9)
        self.assertWithheld(self.assess(short), reason)
        self.assertWithheld(self.assess(replace(short, acquisition_time_s=8.000000004)), reason)

    def test_invalid_extra_recovery_keeps_the_measurement(self):
        independent = tuple(c for c in self.complete.assay.prerequisites if not c.name.endswith(" recovery"))
        for value in (None, -1.0, float("nan")):
            with self.subTest(value=value):
                assessment = self.assess(recovery=replace(RECOVERY, extra_recovery_s=value))
                result = self.assertWithheld(assessment, "extra_recovery_s: Must be finite and nonnegative")
                self.assertEqual((result.I_a, result.I_s, result.R),
                                 (self.complete.result.I_a, self.complete.result.I_s, self.complete.result.R))
                self.assertEqual(assessment.assay.analyte_areas, self.complete.assay.analyte_areas)
                checks = assessment.assay.prerequisites
                self.assertTrue(set(independent).issubset(checks))
                self.assertTrue(all(check.status == "cannot be assessed" for check in checks
                                    if check.name.endswith(" recovery")))

    def test_numeric_processing_error_keeps_its_diagnosis(self):
        with patch("spinivo.pipeline.process", side_effect=FloatingPointError("overflow in spectral transform")):
            self.assertWithheld(self.assess(), "processing: overflow in spectral transform")

    def assess_with_failing_h5(self, function, bounds_index):
        def trap_first_region(*args, **kwargs):
            if args[bounds_index] == self.signals["H5"].bounds_ppm:
                np.divide(1.0, 0.0)
            return DEFAULT

        with np.errstate(divide="raise"), patch(f"spinivo.pipeline.{function.__name__}", wraps=function,
                                                side_effect=trap_first_region):
            return self.assess()

    def test_numeric_region_errors_keep_the_signal_and_diagnosis(self):
        for function, bounds_index, h5_source in ((check_region, 1, None),
                                                  (measure_region, 2, "FID06 H5 analyte region 1")):
            with self.subTest(function=function.__name__):
                assessment = self.assess_with_failing_h5(function, bounds_index)
                result = self.assertWithheld(assessment, "H5 analyte processing: divide by zero encountered in divide")
                self.assertEqual(result.I_s, self.complete.result.I_s)
                self.assertEqual((result.I_a, result.R), (None, None))
                self.assertEqual(assessment.assay.analyte_areas, (None, self.complete.assay.analyte_areas[1]))
                self.assertEqual(assessment.assay.analyte_area_sources, (h5_source, "FID06 I5 analyte region 2"))

    def test_failed_measurement_keeps_its_region_and_pair_checks(self):
        names = ("H5 analyte region 1", "H5/G5 capture_pair")
        region_checks = [check for check in self.complete.assay.prerequisites if check.name in names]
        checks = self.assess_with_failing_h5(measure_region, 2).assay.prerequisites
        self.assertEqual([check for check in checks if check.name in names], region_checks)

    def test_missing_region_keeps_the_declared_region_numbers(self):
        assessment = self.assess(analyte=(replace(self.signals["H5"], bounds_ppm=None), self.signals["I5"]))
        self.assertWithheld(assessment, "H5 analyte region: Missing region")
        self.assertEqual(assessment.assay.analyte_area_sources, (None, "FID06 I5 analyte region 2"))
        names = {check.name for check in assessment.assay.prerequisites}
        self.assertIn("I5 analyte region 2", names)
        self.assertFalse(any(name.startswith("analyte region") for name in names), names)

    def test_numeric_acquisition_errors_mark_only_the_failed_check(self):
        tiny = replace(self.signals["H5"], t1_s=np.float64(1e-320))
        with np.errstate(over="raise"):
            assessment = self.assess(analyte=(tiny, self.signals["I5"]))
        result = self.assertWithheld(assessment, "H5 recovery: tr_s / t1_s must be finite to evaluate recovery")
        self.assertEqual(result.R, self.complete.result.R)
        statuses = {check.name: check.status for check in assessment.assay.prerequisites}
        self.assertEqual((statuses["H5 snr"], statuses["H5 resolution"]), ("pass", "pass"))

    def test_repeated_or_miscounted_signal_names_block_before_measuring(self):
        h5, i5, g5 = (self.signals[name] for name in ("H5", "I5", "G5"))
        for assay, analyte, reason in ((replace(self.assay, analyte_signals=("H5", "H5")), (h5, h5),
                                        "H5 analyte signal: Signal name is repeated in the assay"),
                                       (self.assay, (h5, h5, i5), "H5 analyte signal: Supplied 2 times; declared 1")):
            with self.subTest(reason=reason), patch("spinivo.pipeline.measure_region") as measure:
                self.assertWithheld(self.assess(assay=assay, analyte=analyte, standard=(g5,)), reason)
                measure.assert_not_called()

    def test_missing_sources_block_an_otherwise_passing_assessment(self):
        h5, i5 = (self.signals[name] for name in ("H5", "I5"))
        for changes, reasons in (
            ({"analyte": (replace(h5, t1_source=None), i5)}, ("H5 recovery: t1_source is missing",)),
            ({"analyte": (replace(h5, t1_source="NF"), i5)}, ("H5 recovery: t1_source is missing",)),
            ({"recovery": replace(RECOVERY, calibration_source=None)},
             tuple(f"{name} recovery: calibration_source is missing" for name in CENTERS)),
        ):
            with self.subTest(reasons=reasons):
                self.assertEqual(self.assertWithheld(self.assess(**changes), *reasons).R, self.complete.result.R)

    def test_jeol_fixture_ratio(self):
        references = load_references()
        fixture = references["fixtures"]["jeol"]
        bounds = fixture["regions_ppm"]
        # The FID06 quantities complete the record; only the spectral evidence of the FID02 tube is under test.
        analyte = DeclaredSignal(fixture["numerator"], tuple(bounds[fixture["numerator"]]), None, 2.0, "illustrative")
        standard = DeclaredSignal(fixture["denominator"], tuple(bounds[fixture["denominator"]]), None, 2.0,
                                  "illustrative")
        noise = tuple(references["noise_windows_ppm"]["jeol"])
        path = FIXTURES / fixture["path"]

        assay = replace(self.assay, analyte_signals=(analyte.name,), standard_signals=(standard.name,))

        complete = self.assess(path, assay, (analyte,), (standard,), noise_bounds_ppm=noise)
        self.assertEqual(self.assess(path, assay, (analyte,), (standard,), noise_bounds_ppm=noise), complete)
        result = self.assertWithheld(complete, "DMSO2/MA capture_pair: Missing or unusable analyte capture. "
                                               "Missing or unusable standard capture")
        self.assertAlmostEqual(result.R, fixture["observed_ratio"], delta=1e-9 * fixture["observed_ratio"])
        self.assertFalse(any("recovery:" in reason for reason in result.reasons), result.reasons)


class ExternalPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assay, cls.expected = load_external_case()
        source = FIXTURES / "records" / "illustrative_pair_record.txt"
        cls.assay = replace(cls.assay, pair=replace(cls.assay.pair, source=str(source)))
        cls.raws, cls.signals = {}, {}
        for role in ("sample", "calibrant"):
            run, (name,) = getattr(cls.assay, role), getattr(cls.assay, f"{role}_signals")
            # Matched runs: the same noise-free lines scaled by each fixture area, so R is the area ratio.
            cls.raws[role] = replace(make_raw({"id": run.run_id}, ideal=True), fid=run.area * _fid(((0, 1, 1), *ANCHORS)),
                                     relaxation_delay_s=60.0, pulse_width_us=10.0, pulse_program=run.pulse_sequence,
                                     scans=run.scans, receiver_gain=run.receiver_gain, temperature_k=298.15,
                                     intensity_scale=1.0, storage_scale="JEOL")
            evidence = ComponentEvidence((Line(CARRIER_PPM, run.area, 1.0),), "fixture",
                                         satellite_convention="absent", satellite_source="fixture")
            cls.signals[role] = (DeclaredSignal(name, ppm_bounds((-50, 50)), evidence, 2.0, "T1 log"),)
        cls.complete = cls.assess()

    @classmethod
    def assess(cls, sample_fid=None, calibrant_fid=None, assay=None, **changes):
        inputs = {f"{role}_{name}": value for role in cls.raws for name, value in (
            ("noise_bounds_ppm", ppm_bounds(NOISE_BOUNDS_HZ)), ("recovery", RECOVERY))}
        return assess_external_calibration(sample_fid or cls.raws["sample"], calibrant_fid or cls.raws["calibrant"],
                                           assay or cls.assay, **{**cls.signals, **inputs, **changes})

    def assertBlocked(self, assessment, reason, areas=None):
        result, complete = assessment.result, self.complete.result
        self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
        self.assertIn(reason, result.reasons)
        self.assertEqual((result.J_a, result.J_s, result.R), areas or (complete.J_a, complete.J_s, complete.R))

    def test_matched_runs_pass_with_hand_worked_purity(self):
        assay, result = self.complete.assay, self.complete.result
        self.assertEqual((result.status, result.reasons), ("pass", ()))
        self.assertTrue(all(check.status == "pass" for check in assay.prerequisites))
        self.assertAlmostEqual(result.R, self.expected["R"], delta=1e-9 * self.expected["R"])
        self.assertAlmostEqual(result.value, self.expected["purity"], delta=1e-9)
        self.assertIn("Purity = 97.35 %, u_c = 0.66 %", external_text(external_report(self.complete)))
        with self.assertRaisesRegex(ValueError, "checked pass report"):
            external_report(self.complete.result)
        direct = calculate_external_calibration(assay)
        self.assertEqual((direct.status, direct.reasons, direct.value), ("cannot be assessed", (CALCULATION_ONLY,),
                                                                         result.value))

    def test_repeated_assessments_are_identical(self):
        self.assertEqual(self.assess(), self.complete)
        self.assertEqual(self.assess(assay=self.complete.assay), self.complete)

    def test_baseline_shape_combines_external_tube_terms_in_quadrature(self):
        areas = (BaselineShapeAreas(10, 11, 10), BaselineShapeAreas(20, 20, 16))
        with patch("spinivo.pipeline.baseline_shape_areas", side_effect=areas):
            assessed = self.assess()
        shape = assessed.baseline_shape
        self.assertEqual(shape.status, "computed")
        self.assertEqual(tuple(region.candidate_areas for region in shape.regions), areas)
        self.assertAlmostEqual(shape.relative_standard_uncertainty_percent, 100 * hypot(0.1, 0.2) / sqrt(3))
        self.assertEqual(assessed.result, self.complete.result)
        self.assertEqual(external_report(assessed), external_report(self.complete))

    def test_pair_source_must_be_readable(self):
        pair = self.assay.pair
        missing = replace(self.assay, pair=replace(pair, source=pair.source + ".missing"))
        self.assertBlocked(self.assess(assay=missing), "pair source file is missing")
        source = Path(pair.source)
        dotted = replace(self.assay, pair=replace(pair, source=str(source.parent / "no-such-dir" / ".." / source.name)))
        self.assertBlocked(self.assess(assay=dotted), "pair source file is missing")
        with patch("spinivo.external_calibration.open", side_effect=PermissionError("denied")):
            self.assertBlocked(self.assess(), "pair source file could not be read: denied")
        self.assertEqual(self.assess(), self.complete)

    def test_missing_external_run_is_filled_by_field_name(self):
        absent = replace(self.assay, sample=None)
        with patch.object(ExternalRun, "__init__", autospec=True, side_effect=ExternalRun.__init__) as init:
            assessment = self.assess(assay=absent)
        empty = dict.fromkeys(field.name for field in fields(ExternalRun) if field.init)
        self.assertTrue(all(len(call.args) == 1 for call in init.call_args_list))
        self.assertIn(empty, [call.kwargs for call in init.call_args_list])
        result = assessment.result
        self.assertBlocked(assessment, "sample run: Missing run record", (None, self.complete.result.J_s, None))
        self.assertEqual(result.c_s, self.complete.result.c_s)

    def test_removed_linked_run_drops_its_report_evidence(self):
        linked = self.complete.assay
        original = external_report(calculate_external_calibration(linked))
        for role, other, area in (("sample", "calibrant", "J_s"), ("calibrant", "sample", "J_a")):
            with self.subTest(role=role):
                caller_check = replace(next(check for check in linked.prerequisites if role in check.name))
                assay = replace(linked, prerequisites=(*linked.prerequisites, caller_check), **{role: None})
                report = external_report(calculate_external_calibration(assay))
                self.assertIn(f"{role} run: Missing run record", report["reasons"])
                for field in ("area", "header", "entered"):
                    self.assertIsNone(report[role][field])
                self.assertEqual(report[other], original[other])
                self.assertEqual((report[area], report["c_s"], report["R"]),
                                 (original[area], original["c_s"], None))
                self.assertEqual(report["prerequisites"], [
                    asdict(caller_check), *(check for check in original["prerequisites"] if other in check["name"]),
                    asdict(Prerequisite(f"{role} run", "cannot be assessed", "Missing run record")),
                ])

    def test_pipeline_requires_transfer_and_drift_sources(self):
        for name in ("K", "D"):
            with self.subTest(quantity=name):
                assay = replace(self.assay, **{name: replace(getattr(self.assay, name), source=None)})
                self.assertBlocked(self.assess(assay=assay), f"{name}.source is missing")

    def test_pair_source_expands_only_a_leading_tilde(self):
        root = Path(self.enterContext(TemporaryDirectory()))
        home = root / "home"
        literal, expanded = root / "~" / "pair.txt", home / "pair.txt"
        for path in (literal, expanded):
            path.parent.mkdir()
            path.write_text("Pair evidence", encoding="utf-8")
        sample = replace(self.raws["sample"], source=literal)
        run_id = str(literal)
        assay = replace(self.assay, sample=replace(self.assay.sample, run_id=run_id),
                        pair=replace(self.assay.pair, sample_run_id=run_id, source="./~/pair.txt"))
        with chdir(root), patch.dict(environ, {"HOME": str(home)}):
            self.assertBlocked(self.assess(sample, assay=assay), "pair source names a run file or folder, not a pair record")
            other = replace(assay, pair=replace(assay.pair, source="~/pair.txt"))
            self.assertEqual(self.assess(sample, assay=other).result.status, "pass")

    def test_placeholders_in_nested_evidence_are_missing_before_checks(self):
        signal, = self.signals["sample"]
        components = replace(signal.evidence, target_source="NF", satellite_source="NOT FOUND")
        missing = self.assess(sample=(replace(signal, evidence=components),))
        self.assertRegex("\n".join(missing.result.reasons), r"(?s)Missing target source.*Missing satellite source")
        self.assertEqual(missing.result.status, "cannot be assessed")
        slow = replace(signal, t1_s=1000)
        for changes, field in (
            ({"sample": (replace(slow, t1_source="NF"),)}, "t1_source"),
            ({"sample": (slow,), "sample_recovery": replace(RECOVERY, calibration_source="NF")}, "calibration_source"),
        ):
            with self.subTest(field=field):
                self.assertBlocked(self.assess(**changes), f"sample {signal.name} recovery: {field} is missing")
        assay = replace(self.assay, R_source_ids=("NF",), S=replace(self.assay.S, source_ids=("NOT FOUND",)))
        missing = self.assess(assay=assay)
        for name in ("R", "S"):
            self.assertBlocked(missing, f"{name}.source_ids must contain nonblank source identifiers without outer whitespace")

    def test_header_survives_processing_failure_and_keeps_entered_observations(self):
        sample = replace(self.assay.sample, receiver_gain=31, temperature=None)
        failure = FloatingPointError("overflow in spectral transform")
        with patch("spinivo.pipeline.process", wraps=process, side_effect=[failure, DEFAULT]):
            failed = self.assess(assay=replace(self.assay, sample=sample))
        other = (None, self.complete.result.J_s, None)
        self.assertBlocked(failed, "sample processing: overflow in spectral transform", other)
        self.assertBlocked(failed, "sample receiver_gain: record 31 differs from header 30", other)
        run = external_report(failed)["sample"]
        self.assertEqual((run["scans"], run["receiver_gain"], run["temperature"]), (32, 31, 25.0))
        self.assertEqual((run["header"]["scans"], run["header"]["receiver_gain"]), (32, 30))
        self.assertEqual((run["entered"]["scans"], run["entered"]["receiver_gain"]), (32, 31))
        self.assertNotIn("temperature", run["entered"])
        self.assertEqual(run["header"]["acquisition_time_s"], self.raws["sample"].acquisition_time_s)

    def test_raw_placeholders_are_missing_before_checks_and_header_capture(self):
        raw = replace(self.raws["sample"], pulse_program="NF", storage_scale="NOT FOUND")
        for fid in (raw, raw.source):
            with self.subTest(fid=type(fid)), patch("spinivo.pipeline.read", return_value=raw):
                assessment = self.assess(fid)
                self.assertTrue(any("pulse_program is missing" in reason for reason in assessment.result.reasons))
                header = external_report(assessment)["sample"]["header"]
                self.assertEqual((header["pulse_sequence"], header["storage_scale"]), (None, None))

    def test_report_keeps_entered_scales_separate_from_effective_scales(self):
        sample = replace(self.assay.sample, reader_scale=2, storage_scale="NC=-4")
        linked = self.assess(assay=replace(self.assay, sample=sample))
        run = external_report(self.assess(assay=linked.assay))["sample"]
        self.assertEqual((run["reader_scale"], run["storage_scale"]), (1, "JEOL"))
        self.assertEqual((run["entered"]["reader_scale"], run["entered"]["storage_scale"]), (2, "NC=-4"))

    def test_one_sided_noise_bounds_block_only_the_affected_run_before_measuring(self):
        for role, bounds in (("sample", (1.0,)), ("calibrant", (None, 1.0))):
            with self.subTest(role=role), patch("spinivo.pipeline.measure_region", wraps=measure_region) as measure:
                result = self.assess(**{f"{role}_noise_bounds_ppm": bounds})
                intact = "calibrant" if role == "sample" else "sample"
                self.assertEqual([call.args[1] for call in measure.call_args_list], [self.signals[intact][0].name])
                self.assertEqual(getattr(result.assay, intact), getattr(self.complete.assay, intact))
                self.assertIsNone(getattr(result.assay, role).area)
                self.assertIn(f"{role} noise bounds: Declare two finite ppm bounds", result.result.reasons)

    def test_reference_without_a_source_blocks_only_that_run(self):
        self.assertBlocked(self.assess(sample_reference=(0.5,)), "sample reference: Declare shift_ppm with its source",
                           (None, self.complete.result.J_s, None))

    def test_half_scaled_run_with_sourced_S_keeps_the_purity(self):
        raw = replace(self.raws["sample"], fid=0.5 * self.raws["sample"].fid, receiver_gain=15)
        assay = replace(self.assay, sample=replace(self.assay.sample, receiver_gain=15), S_gains=(15, 30),
                        S=replace(self.assay.S, value=2 * self.assay.S.value))
        scaled = self.assess(raw, assay=assay).result
        self.assertEqual(scaled.status, "pass")
        self.assertAlmostEqual(scaled.value, self.complete.result.value, delta=1e-12)

    def test_missing_evidence_blocks_with_a_named_reason_and_keeps_the_intact_run(self):
        folder = "spectra/sample/10"
        folder_assay = replace(self.assay, sample=replace(self.assay.sample, run_id=folder),
                               pair=replace(self.assay.pair, sample_run_id=folder))
        folder_raw = replace(self.raws["sample"], source=Path(folder))
        self.assertEqual(self.assess(folder_raw, assay=folder_assay).result.status, "pass")
        for changes, reason in (
            ({"assay": replace(self.assay, pair=None)}, "pair record is missing"),
            ({"assay": replace(self.assay, S=replace(self.assay.S, source="NF"))}, "S.source is missing"),
            ({"sample_fid": folder_raw, "assay": replace(folder_assay, pair=replace(folder_assay.pair, source=folder + "/acqus"))},
             "pair source names a run file or folder, not a pair record"),
            ({"calibrant_recovery": replace(RECOVERY, calibrated_flip_angle_deg=30.0)},
             "calibrated_flip_angle_deg differs without K evidence"),
        ):
            with self.subTest(reason=reason):
                self.assertBlocked(self.assess(**changes), reason)
        with patch("spinivo.pipeline.read", side_effect=OSError("corrupt data")):
            unread = self.assess(calibrant_fid=JEOL)
        self.assertBlocked(unread, "calibrant reader: corrupt data",
                           (self.complete.result.J_a, None, None))
        self.assertEqual(unread.assay.sample, self.complete.assay.sample)
        self.assertEqual(external_report(unread)["sample"]["per_scan_area"], self.complete.result.J_a)

    def test_unsupported_input_propagates_for_either_run(self):
        for role in ("sample", "calibrant"):
            with self.subTest(role=role):
                with self.assertRaises(UnsupportedInput) as caught:
                    self.assess(**{f"{role}_fid": JEOL_ARRAY})
                self.assertEqual(caught.exception.reason, "arrayed_or_2d")

    def test_relinking_reads_the_current_header_and_keeps_entered_values(self):
        recovery = RecoveryEvidence(calibrated_90_us=10.0, calibration_source="pulse log")
        first = self.assess(sample_recovery=recovery)
        self.assertEqual(first.result.status, "pass")
        again = self.assess(replace(self.raws["sample"], pulse_width_us=5.0), assay=first.assay, sample_recovery=recovery)
        self.assertBlocked(again, "calibrated_flip_angle_deg differs without K evidence")
        self.assertEqual(again.assay.sample.calibrated_flip_angle_deg, 45.0)
        no_scans = replace(self.raws["sample"], scans=None)
        kept = self.assess(no_scans)
        self.assertBlocked(kept, "sample scans: header has no value; entered 32 is not header evidence")
        run = external_report(kept)["sample"]
        self.assertEqual((run["area"], run["scans"], run["per_scan_area"]),
                         (self.complete.assay.sample.area, 32, self.complete.result.J_a))
        filled = self.assess(assay=replace(self.assay, sample=replace(self.assay.sample, scans=None)))
        self.assertEqual((filled.result.status, filled.assay.sample.scans), ("pass", 32))
        stale = self.assess(no_scans, assay=filled.assay)
        self.assertBlocked(stale, "sample scans is missing", (None, self.complete.result.J_s, None))
        self.assertIsNone(stale.assay.sample.scans)
        entered = filled.entered_assay
        self.assertIsNone(entered.sample.scans)
        reentered = replace(entered, sample=replace(entered.sample, scans=32))
        kept = self.assess(no_scans, assay=reentered)
        self.assertBlocked(kept, "sample scans: header has no value; entered 32 is not header evidence")
        self.assertEqual(kept.entered_assay.sample.scans, 32)


class PipelineInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.internal, cls.external = PipelineTests, ExternalPipelineTests
        cls.internal.setUpClass()
        cls.external.setUpClass()

    def test_raw_source_requires_path_text(self):
        class PlaceholderPath:
            def __fspath__(self):
                return "NF"

        for source in (Path("NF"), False, PlaceholderPath()):
            for method, raw, prefix in ((self.internal, self.internal.raw, ""),
                                        (self.external, self.external.raws["sample"], "sample ")):
                with self.subTest(source=source, method=method.__name__):
                    result = method.assess(replace(raw, source=source)).result
                    self.assertIn(f"{prefix}source: Missing raw source", result.reasons)
                    self.assertEqual((result.status, result.value), ("cannot be assessed", None))

    def test_placeholder_reference_blocks_both_pipelines(self):
        for method, prefix in ((self.internal, ""), (self.external, "sample "), (self.external, "calibrant ")):
            with self.subTest(method=method.__name__, run=prefix):
                key = f"{prefix}reference".replace(" ", "_")
                self.assertEqual(method.assess(**{key: None}), method.complete)
                result = method.assess(**{key: "NF"}).result
                self.assertEqual((result.status, result.value), ("cannot be assessed", None))
                self.assertIn(f"{prefix}reference: Declare shift_ppm with its source", result.reasons)

    def test_dataclass_raw_source_uses_original_path_text(self):
        @dataclass(frozen=True)
        class Source:
            path: str = "NF"

            def __fspath__(self):
                return self.path or "measurement.jdf"

        for method, raw, prefix in ((self.internal, self.internal.raw, ""),
                                    (self.external, self.external.raws["sample"], "sample ")):
            with self.subTest(method=method.__name__):
                result = method.assess(replace(raw, source=Source())).result
                self.assertIn(f"{prefix}source: Missing raw source", result.reasons)
                self.assertEqual((result.status, result.value), ("cannot be assessed", None))

    def test_independent_run_input_errors_are_reported_together(self):
        for method, raw, role in ((self.internal, self.internal.raw, ""),
                                  (self.external, self.external.raws["sample"], "sample")):
            prefix = f"{role} " if role else ""
            key = f"{role}_" if role else ""
            for missing_source, malformed_noise in ((True, False), (False, True), (True, True)):
                with self.subTest(method=method.__name__, missing_source=missing_source, malformed_noise=malformed_noise):
                    fid = replace(raw, source=None) if missing_source else raw
                    changes = {f"{key}reference": "NF"}
                    reasons = [f"{prefix}reference: Declare shift_ppm with its source"]
                    if missing_source:
                        reasons.append(f"{prefix}source: Missing raw source")
                    if malformed_noise:
                        changes[f"{key}noise_bounds_ppm"] = (1.0,)
                        reasons.append(f"{prefix}noise bounds: Declare two finite ppm bounds")
                    result = method.assess(fid, **changes).result
                    self.assertEqual((result.status, result.value), ("cannot be assessed", None))
                    self.assertTrue(set(reasons).issubset(result.reasons), result.reasons)

    def test_missing_raw_source_keeps_acquisition_time_mismatch(self):
        for method, raw, prefix in ((self.internal, self.internal.raw, ""),
                                    (self.external, self.external.raws["sample"], "sample ")):
            with self.subTest(method=method.__name__):
                result = method.assess(replace(raw, source=None, acquisition_time_s=100.0)).result
                self.assertIn(f"{prefix}source: Missing raw source", result.reasons)
                self.assertIn(f"{prefix}acquisition_time_s: 100 s does not match "
                              "complex_points / spectral_width_hz (8 s)", result.reasons)
