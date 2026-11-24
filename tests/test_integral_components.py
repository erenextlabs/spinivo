import unittest
from dataclasses import replace
from math import inf

import numpy as np

from spinivo.integral_checks import PILOT
from spinivo.integral_components import ComponentEvidence, Line, _satellites, check_capture_pair, check_region
from spinivo.processing import process
from tests.integral_fixtures import (
    FOREIGN, NOISE_BOUNDS_HZ, OBSERVE_MHZ, SETTINGS, _fid, component_area, load_cases, make_raw, ppm_bounds,
    target_lines,
)


def _assigned_lines(group):
    return tuple(Line(ppm_bounds((hz,))[0], area, width + SETTINGS.lb_hz) for hz, area, width in group)


class IntegralComponentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        manifest = load_cases()
        cls.cases = {case["id"]: case for case in manifest["cases"] + manifest["variants"]}
        cls.noise = ppm_bounds(NOISE_BOUNDS_HZ)

    def inputs(self, name):
        case = self.cases[name]
        bounds, raw = ppm_bounds(case["bounds"]), make_raw(case)
        spectrum = process(raw, SETTINGS, regions_ppm=(bounds,))

        satellites = case.get("target") == "satellites"
        foreign = () if case.get("withheld") == "extra_line_assignment" else FOREIGN.get(case.get("foreign"), ())
        evidence = ComponentEvidence(
            target_lines=_assigned_lines(target_lines(case)[1:2] if satellites else target_lines(case)),
            target_source="fixture",
            foreign_lines=_assigned_lines(foreign), foreign_source="fixture" if foreign else None,
            satellite_convention="include" if satellites else "absent",
            coupling_hz=130 if satellites else None, satellite_source="fixture",
        )
        return spectrum, bounds, evidence, raw

    def check(self, name, **changes):
        spectrum, bounds, evidence, raw = self.inputs(name)
        return check_region(spectrum, bounds, replace(evidence, **changes), raw=raw, noise_bounds_ppm=self.noise)

    def test_fixture_statuses(self):
        for name, status, cause in (
            ("clean_singlet", "pass", None), ("assigned_triplet", "pass", None), ("satellites_inside", "pass", None),
            ("foreign_line", "flag", "neighbour"), ("foreign_tail", "flag", "neighbour"),
            ("clipped_triplet", "flag", "capture"), ("narrow_singlet", "flag", "capture"),
            ("satellites_outside", "flag", "satellites"),
        ):
            with self.subTest(name=name):
                result = self.check(name)
                self.assertEqual(result.status, status)
                if cause:
                    self.assertIn(cause, dict(result.flags))
                else:
                    self.assertEqual([obs.status for obs in result.observations], ["pass"] * 8)
        result = self.check("unassigned_line")
        self.assertEqual(result.status, "cannot be assessed")
        self.assertIn(("neighbour", "Unexpected peak at 5.07000 ppm"), result.gaps)

    def test_satellite_asymmetry_and_exclusion(self):
        result = self.check("one_satellite_clipped")
        self.assertEqual(result.status, "flag")
        self.assertIn(("satellites", "satellite_asymmetric"), result.flags)
        self.assertLess(result.satellites.value, 0.01)
        self.assertEqual([band.status for band in result.satellites.bands], ["inside", "clipped"])
        for name, status in (("satellites_inside", "flag"), ("satellites_outside", "pass")):
            case = self.cases[name]
            result = self.check(name, satellite_convention="exclude").satellites
            self.assertEqual(result.status, status)
            self.assertAlmostEqual(result.value, component_area(target_lines(case)[::2], case["bounds"]))

    def test_two_parent_satellite_error_exceeds_limit(self):
        bounds_hz = (-60, 60)
        parents = ((-10, 0.5, 1), (10, 0.5, 1))
        satellites = ((-60, 0.00275, 1), (40, 0.00275, 1), (-40, 0.00275, 1), (60, 0.00275, 1))
        evidence = ComponentEvidence(
            target_lines=_assigned_lines(parents), target_source="fixture",
            satellite_convention="include", coupling_hz=100, satellite_source="fixture",
        )
        thresholds = replace(PILOT, max_relative_change=0.001)
        result = _satellites(evidence, OBSERVE_MHZ, ppm_bounds(bounds_hz), (), thresholds)
        expected = 0.011 - component_area(satellites, bounds_hz)
        self.assertAlmostEqual(expected, 0.00286956, places=8)
        self.assertAlmostEqual(result.value, expected, places=14)
        self.assertEqual(result.status, "flag")
        self.assertEqual(result.flags, ("Satellite area error exceeds the limit", "satellite_asymmetric"))
        self.assertEqual([band.status for band in result.bands], ["clipped", "inside", "inside", "clipped"])

    def test_satellite_error_weights_unequal_parents(self):
        bounds_hz = (-70, 60)
        parents = ((0, 1, 1), (10, 3, 3))
        satellites = ((-50, 0.001375, 1), (50, 0.001375, 1), (-40, 0.004125, 3), (60, 0.004125, 3))
        evidence = ComponentEvidence(
            target_lines=_assigned_lines(parents), target_source="fixture", coupling_hz=100, satellite_source="fixture",
        )
        included = component_area(satellites, bounds_hz)
        for convention, expected in (("include", 0.011 - included), ("exclude", included)):
            with self.subTest(convention=convention):
                result = _satellites(replace(evidence, satellite_convention=convention), OBSERVE_MHZ,
                                     ppm_bounds(bounds_hz), (), PILOT)
                self.assertAlmostEqual(result.value, expected, places=14)
                self.assertEqual(result.status, "flag")
                self.assertEqual(result.flags, ("satellite_asymmetric",))
                self.assertEqual([band.status for band in result.bands], ["inside", "inside", "inside", "clipped"])
                for band, (hz, _, width) in zip(result.bands, satellites, strict=True):
                    self.assertAlmostEqual(band.ppm, ppm_bounds((hz,))[0])
                    self.assertEqual(band.fwhm_hz, width + SETTINGS.lb_hz)
                    self.assertAlmostEqual(band.captured, component_area(((hz, 1, width),), bounds_hz), places=13)

    def test_singlet_satellite_inclusion_matches_line_model(self):
        case = self.cases["satellites_inside"]
        result = self.check("satellites_inside").satellites
        self.assertAlmostEqual(result.value, 0.011 - component_area(target_lines(case)[::2], case["bounds"]), places=14)
        self.assertEqual(result.status, "pass")
        self.assertEqual([band.status for band in result.bands], ["inside", "inside"])

    def test_decoupled_satellites_need_no_parent_or_coupling(self):
        evidence = ComponentEvidence(satellite_convention="absent", satellite_source="decoupled control")
        result = _satellites(evidence, OBSERVE_MHZ, ppm_bounds((-60, 60)), (), PILOT)
        self.assertEqual((result.status, result.value, result.bands), ("pass", 0.0, ()))

    def test_clean_assigned_triplet_with_satellites_passes(self):
        case = dict(self.cases["assigned_triplet"], snr=30000)
        parents = target_lines(case)
        satellites = tuple((hz + offset, area * 0.0055, width)
                           for hz, area, width in parents for offset in (-65, 65))
        raw = make_raw(case)
        raw = replace(raw, fid=raw.fid + _fid(satellites))
        bounds = ppm_bounds(case["bounds"])
        spectrum = process(raw, SETTINGS, regions_ppm=(bounds,))
        evidence = ComponentEvidence(
            target_lines=_assigned_lines(parents), target_source="fixture",
            satellite_convention="include", coupling_hz=130, satellite_source="fixture",
        )
        result = check_region(spectrum, bounds, evidence, raw=raw, noise_bounds_ppm=self.noise)
        self.assertEqual(result.status, "pass", result.flags + result.gaps)
        self.assertEqual(result.neighbour.status, "pass")
        self.assertAlmostEqual(result.satellites.value, 0.011 - component_area(satellites, case["bounds"]), places=14)
        self.assertEqual([band.status for band in result.satellites.bands], ["inside"] * 6)

    def test_missing_evidence_and_flag_precedence(self):
        spectrum, bounds, evidence, raw = self.inputs("foreign_line")
        result = check_region(spectrum, bounds, evidence, raw=raw)
        self.assertEqual(result.status, "flag")
        self.assertIn("neighbour", dict(result.flags))
        self.assertIn(("neighbour", "Missing noise window"), result.gaps)
        self.assertIn("Missing foreign source", self.check("foreign_line", foreign_source=" ").neighbour.reason)
        spectrum, bounds, evidence, raw = self.inputs("clean_singlet")
        for changed, noise, gap in (
            (replace(evidence, target_source=" \t"), self.noise, ("capture", "Missing target source")),
            (replace(evidence, satellite_source=""), self.noise, ("satellites", "Missing satellite source")),
            (replace(evidence, satellite_convention=None), self.noise, ("satellites", "Missing satellite convention")),
            (evidence, None, ("neighbour", "Missing noise window")),
        ):
            result = check_region(spectrum, bounds, changed, raw=raw, noise_bounds_ppm=noise)
            self.assertEqual(result.status, "cannot be assessed")
            self.assertIn(gap, result.gaps)
        self.assertEqual(check_region(spectrum, bounds, evidence, raw=raw, noise_bounds_ppm=self.noise).status, "pass")

    def test_pair_capture(self):
        clean = self.check("clean_singlet")
        self.assertEqual(check_capture_pair(clean, clean).status, "pass")
        self.assertEqual(check_capture_pair(self.check("narrow_singlet"), clean).status, "flag")
        total_loss = replace(clean, capture=replace(clean.capture, value=1.0))
        flagged = check_capture_pair(total_loss, clean)
        self.assertEqual((flagged.status, flagged.value), ("flag", 1.0))
        for missing in (None, 1.0, 1.1):
            result = check_capture_pair(clean, replace(clean, capture=replace(clean.capture, value=missing)))
            self.assertEqual(result.status, "cannot be assessed")
            self.assertIn("standard capture", result.reason)
        limit = replace(clean, capture=replace(clean.capture, value=PILOT.max_pair_capture_mismatch))
        full = replace(clean, capture=replace(clean.capture, value=0.0))
        self.assertEqual(check_capture_pair(limit, full).status, "pass")

    def test_unexpected_peak_limit_equality(self):
        spectrum, bounds, evidence, raw = self.inputs("clean_singlet")
        noise = (self.noise[0] <= spectrum.ppm) & (spectrum.ppm <= self.noise[1])
        data = np.where(noise, spectrum.data, 0)
        cutoff = PILOT.unexpected_peak_sd * np.std(data.real[noise], ddof=1)
        peak = np.argmin(abs(spectrum.hz - 80))
        for height, status in ((cutoff, "pass"), (np.nextafter(cutoff, inf), "cannot be assessed")):
            data[peak] = height
            result = check_region(replace(spectrum, data=data), bounds, evidence, raw=raw, noise_bounds_ppm=self.noise)
            self.assertEqual(result.neighbour.status, status)
