import unittest
from dataclasses import replace
from math import isfinite

import numpy as np

from spinivo import read
from spinivo.processing import QUANTITATIVE_SETTINGS, ProcessingSettings, integrate, process
from tests.test_processing import BRUKER_REFERENCE_SHIFT_PPM, JEOL_REFERENCE_SHIFT_PPM
from tests.test_readers import FIXTURES
from tools.compare_references import load_references, measure_regions


class ProcessingReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        references = load_references()
        cls.references = references
        cls.tolerance = references["relative_tolerance"]
        cls.observed_ratios = {}
        cls.residuals = {}
        for vendor, shift, names in (
            ("jeol", JEOL_REFERENCE_SHIFT_PPM, ("vendor", "workbook")),
            ("bruker", BRUKER_REFERENCE_SHIFT_PPM, ("vendor",)),
        ):
            fixture = references["fixtures"][vendor]
            regions = measure_regions(
                FIXTURES / fixture["path"], fixture["regions_ppm"], references,
                processing_override=fixture.get("processing"),
                reference_shift_ppm=shift, reference_source=fixture["reference_source"],
            )
            for region in regions.values():
                if region.area is None or not isfinite(region.area) or region.area <= 0:
                    raise AssertionError(f"{region.source}: invalid area ({region.area_reason or region.area})")
            ratio = regions[fixture["numerator"]].area / regions[fixture["denominator"]].area
            cls.observed_ratios[vendor] = (ratio, fixture["observed_ratio"])
            if vendor == "bruker":
                top, bottom = (regions[fixture[key]].baseline_region for key in ("numerator", "denominator"))
                shifted = top.half_block_area / bottom.half_block_area
                cls.bruker_sensitivity = shifted, (shifted - ratio) / ratio, fixture["baseline_sensitivity"]
            for name in names:
                reference = fixture["references"][name]
                expected = reference["ratio"]
                residual = (ratio - expected) / expected
                message = f"{vendor} {name}: reference={expected:.8g}, ours={ratio:.8g}, residual={residual:+.4%}"
                print(message, flush=True)
                cls.residuals[vendor, name] = (residual, reference["relative_residual"], message)

    def assert_ratio(self, vendor, reference):
        residual, _, message = self.residuals[vendor, reference]
        self.assertLessEqual(abs(residual), self.tolerance, message)

    def test_recorded_measurements(self):
        for vendor, (ratio, recorded) in self.observed_ratios.items():
            with self.subTest(vendor=vendor, field="observed_ratio"):
                self.assertAlmostEqual(ratio, recorded, delta=1e-9 * abs(recorded))
        for (vendor, name), (residual, recorded, message) in self.residuals.items():
            with self.subTest(vendor=vendor, reference=name, field="relative_residual"):
                self.assertAlmostEqual(residual, recorded, delta=1e-9 * abs(recorded), msg=message)

    def test_quantitative_settings(self):
        self.assertEqual(ProcessingSettings(**self.references["processing"]), QUANTITATIVE_SETTINGS)

    def test_rubric_settings(self):
        settings = ProcessingSettings(**self.references["rubric_processing"])
        self.assertEqual(settings.lb_hz, 1)
        self.assertFalse(settings.baseline_bridge)
        self.assertFalse(settings.lock_spike_anchors)

    def test_bruker_processing_override(self):
        override = self.references["fixtures"]["bruker"]["processing"]
        self.assertEqual(override["lb_hz"], 1)
        self.assertIsNone(override["phase_lb_hz"])

    def test_jeol_vendor_ratio(self):
        self.assert_ratio("jeol", "vendor")

    def test_jeol_workbook_ratio(self):
        self.assert_ratio("jeol", "workbook")

    def test_bruker_vendor_ratio(self):
        self.assert_ratio("bruker", "vendor")

    def test_bruker_half_block_sensitivity(self):
        shifted, change, recorded = self.bruker_sensitivity
        self.assertAlmostEqual(shifted, recorded["half_block_ratio"], delta=1e-9 * shifted)
        self.assertAlmostEqual(change, recorded["relative_ratio_change"], delta=1e-9 * abs(change))
        self.assertGreater(abs(change), self.tolerance)
        self.assertEqual(recorded["assessment"], "baseline sensitive")


def _vendor_region_scale(raw, scale):
    hz_per_ppm = raw.observe_frequency_mhz
    if raw.vendor == "jeol":
        factor = scale["delta_zero_fill_factor"] * raw.complex_points
    else:
        nc = int(raw.storage_scale.removeprefix("NC="))
        factor = 2 ** (scale["NC_proc"] - nc)
    return hz_per_ppm / raw.spectral_width_hz * factor


def _bruker_vendor_phase(raw, procs):
    group_delay = raw.digital_filter["GRPDLY"]
    # Convert TopSpin's left-edge phase to Spinivo's carrier pivot after delay removal.
    left = (-procs["PHC0"] - 180 * group_delay + 180) % 360 - 180
    p1 = -procs["PHC1"]
    return left + p1 / 2, p1


class VendorRegionAreaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixtures = load_references()["fixtures"]
        cls.raw = {vendor: read(FIXTURES / fixture["path"]) for vendor, fixture in cls.fixtures.items()}

    def region_errors(self, vendor, *, autophase=False, gain=1.0):
        fixture = self.fixtures[vendor]
        reference = fixture["vendor_region_areas"]
        raw = self.raw[vendor]
        shift = JEOL_REFERENCE_SHIFT_PPM if vendor == "jeol" else BRUKER_REFERENCE_SHIFT_PPM
        settings = replace(
            QUANTITATIVE_SETTINGS, smooth_baseline=False, reference_shift_ppm=shift,
            reference_source=fixture["reference_source"], **fixture.get("processing", {}),
        )
        if vendor == "bruker" and not autophase:
            p0, p1 = _bruker_vendor_phase(raw, reference["procs"])
            settings = replace(settings, zero_fill_factor=2, autophase_method=None, baseline_windows_ppm=(),
                               p0_deg=p0, p1_deg=p1, pivot_ppm=None)
        scale = _vendor_region_scale(raw, reference["scale"])
        if gain != 1.0:
            raw = replace(raw, fid=raw.fid * gain)
        regions = reference["regions"]
        windows = tuple(tuple(region["bounds_ppm"]) for region in regions.values())
        spectrum = process(raw, settings, regions_ppm=windows)
        for (name, region), bounds in zip(regions.items(), windows):
            integral = integrate(spectrum, bounds)
            ends = np.interp(bounds, spectrum.ppm[::-1], spectrum.data.real[::-1])
            lo_hz, hi_hz = integral.bounds_hz
            area = integral.area - (hi_hz - lo_hz) * ends.mean()
            error = area / (scale * region["vendor_area"]) - 1
            mode = settings.autophase_method or "vendor phase"
            message = f"{vendor} {name} ({mode}, FID gain {gain:.2f}): error={error:+.6%}"
            print(message, flush=True)
            yield name, region, error, message

    def assert_region_areas(self, vendor, *, autophase=False):
        tolerance = self.fixtures[vendor]["vendor_region_areas"]["rtol"]
        for name, region, error, message in self.region_errors(vendor, autophase=autophase):
            with self.subTest(region=name):
                self.assertLessEqual(abs(error), tolerance, message)
                if not autophase:
                    recorded = region["relative_error"]
                    self.assertAlmostEqual(error, recorded, delta=1e-9 * abs(recorded), msg=message)

    def test_jeol_vendor_region_areas(self):
        self.assert_region_areas("jeol")

    def test_bruker_vendor_region_areas(self):
        self.assert_region_areas("bruker")

    def test_bruker_vendor_region_areas_autophase(self):
        self.assert_region_areas("bruker", autophase=True)

    def test_vendor_region_areas_detect_gain_error(self):
        for vendor in ("jeol", "bruker"):
            tolerance = self.fixtures[vendor]["vendor_region_areas"]["rtol"]
            for name, _, error, message in self.region_errors(vendor, gain=1.02):
                with self.subTest(vendor=vendor, region=name):
                    self.assertTrue(isfinite(error), message)
                    self.assertGreater(abs(error), tolerance, message)

    def test_vendor_region_scale(self):
        for vendor, raw in self.raw.items():
            scale = self.fixtures[vendor]["vendor_region_areas"]["scale"]
            expected = scale["expected"]
            with self.subTest(vendor=vendor):
                self.assertAlmostEqual(_vendor_region_scale(raw, scale), expected, delta=1e-12 * abs(expected))
        reference = self.fixtures["bruker"]["vendor_region_areas"]
        phase = _bruker_vendor_phase(self.raw["bruker"], reference["procs"])
        for name, actual in zip(("p0_deg", "p1_deg"), phase):
            expected = reference["expected_phase"][name]
            with self.subTest(phase=name):
                self.assertAlmostEqual(actual, expected, delta=1e-12 * abs(expected))
