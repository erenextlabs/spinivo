import unittest
from math import isfinite

from spinivo.processing import QUANTITATIVE_SETTINGS, ProcessingSettings
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
