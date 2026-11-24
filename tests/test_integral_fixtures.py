import csv
import os
import subprocess
import sys
import unittest
from collections import Counter, defaultdict
from contextlib import redirect_stdout
from dataclasses import replace
from io import StringIO
from math import atan, pi
from pathlib import Path
from unittest.mock import patch

import numpy as np

from spinivo import read
from spinivo.integral_checks import check_baseline_recovery
from spinivo.integral_components import check_region
from spinivo.processing import QUANTITATIVE_SETTINGS, integrate, process
from spinivo.scoring import Prediction, score_integrals
from tests import test_scoring
from tests.integral_fixtures import (
    DEFAULT_BOUNDS_HZ, FIXTURES, FOREIGN, SEEDS, SETTINGS, component_area,
    declared_evidence, load_cases, make_raw, ppm_bounds, processing_regions, target_lines,
)
from tools.compare_references import DEPOSIT, EXTRACTED, load_references, measure_regions, report_checks

PROJECT = FIXTURES.parents[1]


def quantitative_target_area(case):
    low, high = case["bounds"]
    return sum(area * (atan(2 * (high - center) / width) - atan(2 * (low - center) / width)) / pi
               for center, area, width in target_lines(case))


class IntegralFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = load_cases()
        cls.cases = {case["id"]: case for case in cls.manifest["cases"]}

    def test_manifest(self):
        self.assertEqual(len(self.cases), len(self.manifest["cases"]))
        self.assertEqual(Counter((case["kind"], case["label"]) for case in self.cases.values()), {
            ("synthetic", "reliable"): 3, ("synthetic", "unreliable"): 8, ("synthetic", "indeterminate"): 2,
            ("real", "unreliable"): 3, ("real", "indeterminate"): 2,
        })
        for name, case in self.cases.items():
            self.assertTrue(case["source"] and case["evidence"], name)
            self.assertEqual(bool(case["causes"]), case["label"] != "reliable", name)
        candidates = self.manifest["candidates"]
        self.assertEqual([row["path"] for row in candidates], [f"RCE/Accuracy_Sample0{i}" for i in range(1, 4)])
        self.assertTrue(all(row["scored_regions"] == 0 for row in candidates))

    def test_synthetic_processed_areas(self):
        omitted = {"clipped_triplet", "satellites_outside", "narrow_singlet"}
        for name, case in self.cases.items():
            if case["kind"] != "synthetic" or name == "low_snr":
                continue
            with self.subTest(case=name):
                regions = processing_regions(case)
                seeds = SEEDS[:1] if case["label"] == "unreliable" else SEEDS
                spectra = [process(make_raw(case, seed), SETTINGS, regions_ppm=regions) for seed in seeds]
                self.assert_edge_trials(spectra[0], case)
                areas = np.array([integrate(s, regions[0]).area for s in spectra])
                intended = DEFAULT_BOUNDS_HZ if name in omitted else case["bounds"]
                errors = areas / component_area(target_lines(case), intended) - 1
                if case["label"] != "unreliable":
                    self.assertGreater(areas.mean(), 0)
                    self.assertLessEqual(areas.std(ddof=1) / areas.mean(), 0.01)
                    changes = [s.smooth_baseline_fit.regions[0].half_block_relative_change for s in spectra]
                    self.assertLessEqual(np.max(np.abs(np.array(changes, dtype=float))), 0.01)
                if case["label"] == "reliable":
                    self.assertLessEqual(np.max(np.abs(errors)), 0.005)
                    self.assertLessEqual(areas.std(ddof=1) / areas.mean(), 0.003)
                    self.assertTrue(all(s.smooth_baseline_fit.selection.branch == "constant" for s in spectra))
                elif case["label"] == "unreliable":
                    signed_error = -errors[0] if name in omitted or name == "local_phase" else errors[0]
                    self.assertGreater(signed_error, 0.01)

    def test_tail_and_local_corruptions(self):
        tail = self.cases["foreign_tail"]
        self.assertGreater(FOREIGN[tail["foreign"]][0][0], tail["bounds"][1])
        truth = component_area(target_lines(tail), tail["bounds"])
        self.assertGreater(component_area(FOREIGN[tail["foreign"]], tail["bounds"]) / truth, 0.02)

        phase = self.cases["local_phase"]
        regions = processing_regions(phase)
        control = process(make_raw(phase, ideal=True), SETTINGS, regions_ppm=regions)
        rotated = process(make_raw(phase), SETTINGS, regions_ppm=regions)
        truth = component_area(target_lines(phase), phase["bounds"])
        self.assertAlmostEqual(integrate(control, regions[0]).area / truth, 1, delta=0.005)
        self.assertAlmostEqual(phase["phase_deg"] + rotated.settings.p0_deg, 20, delta=0.5)

        baseline = self.cases["local_baseline"]
        regions = processing_regions(baseline)
        clean = process(make_raw(self.cases["clean_singlet"]), SETTINGS, regions_ppm=regions)
        raised = process(make_raw(baseline), SETTINGS, regions_ppm=regions)
        truth = component_area(target_lines(baseline), baseline["bounds"])
        extra = (integrate(raised, regions[0]).area - integrate(clean, regions[0]).area) / truth
        self.assertAlmostEqual(extra, 0.03, delta=0.005)

    def test_low_snr_repeatability(self):
        case = self.cases["low_snr"]
        regions = processing_regions(case)
        control_settings = replace(SETTINGS, smooth_baseline=False)
        areas, controls = [], []
        for seed in SEEDS:
            raw = make_raw(case, seed)
            spectrum = process(raw, SETTINGS, regions_ppm=regions)
            control = process(raw, control_settings, regions_ppm=regions)
            areas.append(integrate(spectrum, regions[0]).area)
            controls.append(integrate(control, regions[0]).area)
            if seed == SEEDS[0]:
                self.assert_edge_trials(spectrum, case)
        self.assertGreater(min(areas + controls), 0)
        self.assertGreater(np.std(controls, ddof=1) / np.mean(controls), 0.01)
        cv = np.std(areas, ddof=1) / np.mean(areas)
        self.assertGreater(cv, 0.01)
        self.assertAlmostEqual(cv, case["area_cv"], delta=1e-7)

    def test_indeterminate_cases_withhold_one_item(self):
        for name, control_name, withheld in (("no_noise_window", "clean_singlet", "noise_window_hz"),
                                              ("unassigned_line", "foreign_line", "extra_line_assignment")):
            case, control = self.cases[name], self.cases[control_name]
            self.assertEqual(case["bounds"], control["bounds"], name)
            np.testing.assert_array_equal(make_raw(case).fid, make_raw(control).fid, err_msg=name)
            evidence = declared_evidence(control)
            del evidence[withheld]
            self.assertEqual(declared_evidence(case), evidence, name)

    def test_satellite_bounds(self):
        variant, = self.manifest["variants"]
        self.assertEqual((variant["label"], variant["causes"], variant["checker_status"], variant["purity_eligible"]),
                         ("unreliable", ["satellites"], "flagged", False))
        case = self.cases[variant["parent"]]
        raw = make_raw(variant)
        np.testing.assert_array_equal(raw.fid, make_raw(case).fid)
        left, _, right = target_lines(case)
        self.assertLess(variant["bounds"][0], left[0])
        self.assertEqual(variant["bounds"][1], right[0])
        satellites = (left, right)
        full = component_area(satellites, case["bounds"])
        both_outside = self.cases["satellites_outside"]["bounds"]
        self.assertGreater(full - component_area(satellites, both_outside), 0.01)
        self.assertAlmostEqual(full - component_area(satellites, variant["bounds"]), 0.00275, delta=0.0001)
        # Evaluate the changed edge on the nominal trace with processing held fixed.
        spectrum = process(raw, SETTINGS, regions_ppm=processing_regions(case))
        self.assert_edge_trials(spectrum, variant)
        bounds = ppm_bounds(variant["bounds"])
        truth = component_area(target_lines(variant), variant["bounds"])
        self.assertAlmostEqual(integrate(spectrum, bounds).area / truth, 1, delta=0.005)

    def assert_edge_trials(self, spectrum, case):
        low, high = ppm_bounds(case["bounds"]) if case["kind"] == "synthetic" else case["bounds"]
        step = max(spectrum.spectral_width_hz / spectrum.acquired_points / spectrum.observe_frequency_mhz,
                   0.05 * (high - low))
        area = integrate(spectrum, (low, high)).area
        self.assertGreater(area, 0)
        trials = ((low + step, high), (low - step, high), (low, high - step), (low, high + step))
        changes = [integrate(spectrum, trial).area / area - 1 for trial in trials]
        if case["kind"] == "real":
            recorded = [case["edge_step_ppm"], *case["edge_area_relative_changes"]]
            np.testing.assert_allclose([step, *changes], recorded, atol=5e-9, rtol=0, err_msg=case["id"])
        self.assertEqual("bounds" in case["causes"], max(map(abs, changes)) > 0.01, case["id"])
        if "bounds" in case["causes"]:
            self.assertEqual(case["label"], "unreliable")

    def test_committed_real_windows(self):
        groups = defaultdict(list)
        for case in self.cases.values():
            if case["kind"] == "real" and not case.get("external"):
                groups[case["source"]].append(case)
        self.assertEqual(sorted(map(len, groups.values())), [2, 2])
        for source, cases in groups.items():
            references = {(case.get("reference_shift_ppm", 0), case.get("reference_source")) for case in cases}
            (shift, provenance), = references
            settings = replace(SETTINGS, reference_shift_ppm=shift, reference_source=provenance)
            regions = tuple(tuple(case["bounds"]) for case in cases)
            spectrum = process(read(PROJECT / source), settings, regions_ppm=regions)
            self.assertEqual(spectrum.smooth_baseline_fit.selection.branch, "smooth")
            for case, region in zip(cases, spectrum.smooth_baseline_fit.regions):
                self.assert_edge_trials(spectrum, case)
                if "half_block_relative_change" in case:
                    self.assertAlmostEqual(region.half_block_relative_change,
                                           case["half_block_relative_change"], delta=5e-7)
                expected = {
                    "bruker_1_15": (1.0396, 0.00003521, 0.03164098, "baseline_sensitive"),
                    "bruker_7_64": (0.9791, 0.00162963, 0.00009996, None),
                }.get(case["id"])
                if expected is not None:
                    midpoint, position, placement, flag = expected
                    recovery = check_baseline_recovery(spectrum, region.bounds_ppm)
                    self.assertAlmostEqual(recovery.value, midpoint, delta=0.0002)
                    self.assertAlmostEqual(recovery.position_sensitivity, position, delta=5e-7)
                    self.assertAlmostEqual(recovery.placement_sensitivity, placement, delta=5e-7)
                    self.assertEqual(recovery.flag, flag)
                    self.assertIsNone(recovery.reason)

    def test_external_window_inventory(self):
        case = self.cases["fid05_widened_dmso2"]
        with (PROJECT / "data/manifest.csv").open(encoding="utf-8", newline="") as inventory:
            members = {row["path"] for row in csv.DictReader(inventory) if row["doi"] == "10.7910/DVN/QFY6L9"}
        self.assertTrue({case["source"], case["source"].replace("-1-1.jdf", "-1-2.jdf")} <= members)

    def test_external_real_window(self):
        root = os.environ.get("SPINIVO_NISHIZAKI_ROOT")
        if not root:
            self.skipTest("SPINIVO_NISHIZAKI_ROOT is unset")
        case = self.cases["fid05_widened_dmso2"]
        bounds = tuple(case["bounds"])
        spectrum = process(read(Path(root) / case["source"]), SETTINGS, regions_ppm=(bounds,))
        self.assert_edge_trials(spectrum, case)
        target = integrate(spectrum, (2.93, 3.04)).area
        self.assertGreater(target, 0)
        self.assertGreater(integrate(spectrum, (3.20, bounds[1])).area / target, 2)

    def test_external_pilot_ratios(self):
        cached = PROJECT / ".cache/veri/nishizaki_QFY6L9" / EXTRACTED.parent
        root = Path(os.environ.get("SPINIVO_NISHIZAKI_ROOT") or cached)
        if root.name != DEPOSIT:
            root /= DEPOSIT
        if not root.is_dir():
            self.skipTest("Extracted Nishizaki pilot data is absent")

        references = load_references()
        limit = references["relative_tolerance"]
        self.assertEqual(limit, 0.01)
        residuals = {}
        for pilot in references["pilots"]:
            measured_areas, reference_areas = {}, {}
            for run in pilot["runs"]:
                regions = run["regions"]
                windows = {cell: region["bounds_ppm"] for cell, region in regions.items()}
                measured = measure_regions(root / run["path"], windows, references)
                for cell, result in measured.items():
                    area = result.area
                    self.assertTrue(area is not None and np.isfinite(area) and area > 0, result)
                    measured_areas[cell] = area
                    reference_areas[cell] = regions[cell]["area"]
            for top, bottom in pilot["ratios"]:
                measured_ratio = measured_areas[top] / measured_areas[bottom]
                expected_ratio = reference_areas[top] / reference_areas[bottom]
                residuals[top, bottom] = measured_ratio / expected_ratio - 1
        self.assertEqual(len(residuals), 7)
        failures = {pair for pair, error in residuals.items() if abs(error) > limit}
        self.assertEqual(failures, {("I149", "H146"), ("H5", "F4"), ("I5", "F4")}, residuals)


class QuantitativeIntegralFixtureTests(unittest.TestCase):
    def test_synthetic_scores(self):
        manifest = load_cases()
        cases = {case["id"]: case for case in manifest["cases"] if case["kind"] == "synthetic"}
        # Use the same evidence and effective linewidth recipe as the historical scoring tests.
        with patch.object(test_scoring, "SETTINGS", QUANTITATIVE_SETTINGS):
            prepared, version = test_scoring.prepare_fixture([*cases.values(), *manifest["variants"]])
        results = {name: check_region(row["spectrum"], row["bounds"], row["components"],
                                      raw=row["raw"], noise_bounds_ppm=row["noise_bounds_ppm"])
                   for name, row in prepared.items()}
        self.assertEqual({result.thresholds.version for result in results.values()}, {"pilot-3"})
        predictions = {name: Prediction.from_region(result) for name, result in results.items()}
        cohorts = score_integrals(manifest, predictions, processing_version=version)["cohorts"]
        synthetic = cohorts["synthetic"]
        self.assertEqual((synthetic["true_positive"], synthetic["false_positive"], synthetic["false_negative"]),
                         (7, 0, 1))
        self.assertEqual((synthetic["precision"], synthetic["recall"], synthetic["assessed"], synthetic["labelled"]),
                         (7 / 7, 7 / 8, 11, 13))
        self.assertEqual(cohorts["variant"]["true_positive"], 1)

        raised, clean = (results[name] for name in ("local_baseline", "clean_singlet"))
        self.assertEqual(raised.status, "pass")
        extra = (raised.nominal_area - clean.nominal_area) / quantitative_target_area(cases["local_baseline"])
        self.assertAlmostEqual(extra, 0.03007, delta=5e-6)
        narrow = results["narrow_singlet"]
        self.assertEqual((narrow.status, predictions["narrow_singlet"].flags), ("flag", ("capture",)))
        self.assertAlmostEqual(narrow.capture.value, 0.15596, delta=5e-6)
        edge = max((trial.relative_change for trial in narrow.bounds.details), key=abs)
        self.assertAlmostEqual(edge, -0.00981, delta=5e-6)

    def test_reliable_controls(self):
        for case in load_cases()["cases"]:
            if case["kind"] != "synthetic" or case["label"] != "reliable":
                continue
            with self.subTest(case=case["id"]):
                regions = processing_regions(case)
                areas = []
                for seed in SEEDS:
                    spectrum = process(make_raw(case, seed), QUANTITATIVE_SETTINGS, regions_ppm=regions)
                    areas.append(integrate(spectrum, regions[0]).area)
                    half_block = spectrum.smooth_baseline_fit.regions[0].half_block_relative_change
                    self.assertTrue(half_block is not None and np.isfinite(half_block), seed)
                    self.assertLessEqual(abs(half_block), 0.01, seed)
                self.assertGreater(min(areas), 0)
                for truth in (quantitative_target_area(case), component_area(target_lines(case), case["bounds"])):
                    self.assertLessEqual(np.max(np.abs(np.array(areas) / truth - 1)), 0.005)
                self.assertLessEqual(np.std(areas, ddof=1) / np.mean(areas), 0.003)

    def test_committed_real_windows_cli(self):
        completed = subprocess.run(
            [sys.executable, str(PROJECT / "tools/compare_references.py"), "--checks", "--quantitative", str(FIXTURES)],
            cwd=PROJECT, check=True, capture_output=True, text=True,
        )
        reports = completed.stdout.strip().split("\nlabel=")
        expected = (
            ("fid02_dmso2", "CANNOT BE ASSESSED", "none", "indeterminate"),
            ("fid02_ma", "CANNOT BE ASSESSED", "none", "indeterminate"),
            ("bruker_1_15", "FLAG", "baseline,bounds,phase", "unreliable"),
            ("bruker_7_64", "FLAG", "phase", "unreliable"),
        )
        self.assertEqual(len(reports), len(expected) + 1)
        for i, (name, status, causes, label) in enumerate(expected):
            with self.subTest(case=name):
                heading = next(line for line in reports[i].splitlines() if line.startswith(f"{name} "))
                self.assertTrue(heading.endswith(f"ppm: {status}"), heading)
                self.assertEqual(reports[i + 1].splitlines()[0], label)
                self.assertRegex(reports[i], rf"(?m)^  causes={causes} threshold_version=pilot-3 "
                                 r"midpoint_recovery=\d+\.\d+% P=\d+\.\d+pp C=\d+\.\d+% "
                                 r"half_block_change=[+-]\d+\.\d+%$")
                self.assertRegex(reports[i], r"(?m)^  edge_moves: low_out=[+-]\d+\.\d+% low_in=[+-]\d+\.\d+% "
                                 r"high_in=[+-]\d+\.\d+% high_out=[+-]\d+\.\d+%$")

    def test_external_real_window(self):
        cached = PROJECT / ".cache/veri/nishizaki_QFY6L9" / EXTRACTED.parent
        root = Path(os.environ.get("SPINIVO_NISHIZAKI_ROOT") or cached)
        if root.name != DEPOSIT:
            root /= DEPOSIT
        case = next(case for case in load_cases()["cases"] if case["id"] == "fid05_widened_dmso2")
        if not (root / Path(case["source"]).relative_to(DEPOSIT)).is_file():
            self.skipTest("Extracted Nishizaki FID05 data is absent")
        output = StringIO()
        with redirect_stdout(output):
            report_checks(root, load_references(), quantitative=True)
        report = output.getvalue().split("\nfid05_widened_dmso2 ")[1]
        self.assertTrue(report.splitlines()[0].endswith("ppm: CANNOT BE ASSESSED"))
        self.assertIn("causes=none threshold_version=pilot-3", report)
        self.assertTrue(report.endswith("label=unreliable\n"))
