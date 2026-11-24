import unittest
from dataclasses import replace
from math import inf, nan
from unittest.mock import patch

import numpy as np

from spinivo import read
from spinivo.integral_checks import check_region_stability
from spinivo.measurement import measure_region
from spinivo.processing import ProcessingSettings, integrate, process
from tests.test_readers import JEOL


REGIONS = ((2.90, 3.01), (6.17, 6.28))
NOISE = (1.55, 1.85)
SETTINGS = ProcessingSettings(autophase_method="lineshape", smooth_baseline=True, pivot_ppm=4.25)


class PhaseTransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = read(JEOL)

    def process(self, settings):
        return process(self.raw, settings, regions_ppm=REGIONS)

    def two_pass(self, settings):
        phase = self.process(replace(settings, lb_hz=settings.phase_lb_hz, phase_lb_hz=None))
        manual = replace(phase.settings, lb_hz=settings.lb_hz, autophase_method=None)
        return phase, self.process(manual)

    def assert_trace_matches(self, transferred, expected):
        np.testing.assert_allclose(transferred.data, expected.data, rtol=1e-10, atol=0)
        np.testing.assert_allclose(transferred.smooth_baseline_curve, expected.smooth_baseline_curve,
                                   rtol=1e-10, atol=0)
        np.testing.assert_allclose([integrate(transferred, bounds).area for bounds in REGIONS],
                                   [integrate(expected, bounds).area for bounds in REGIONS], rtol=1e-10, atol=0)
        self.assertEqual(transferred.smooth_baseline_fit, expected.smooth_baseline_fit)

    def test_matches_two_pass(self):
        for method, lb, phase_lb in (("lineshape", 0, 1), ("acme", 0, 1), ("acme", 1, 0)):
            with self.subTest(method=method, lb_hz=lb, phase_lb_hz=phase_lb):
                settings = replace(SETTINGS, autophase_method=method, lb_hz=lb, phase_lb_hz=phase_lb)
                phase, expected = self.two_pass(settings)
                transferred = self.process(settings)
                self.assert_trace_matches(transferred, expected)
                self.assertEqual(transferred.settings,
                                 replace(expected.settings, autophase_method=method, phase_lb_hz=phase_lb))
                self.assertEqual(transferred.lineshape_phase_fit, phase.lineshape_phase_fit)

    def test_acme_uses_broadened_peak_for_polarity(self):
        settings = replace(SETTINGS, autophase_method="acme", phase_lb_hz=1)
        # At -102 degrees the 1 Hz peak is negative and the 0 Hz peak is positive.
        with patch("spinivo.processing.autops", return_value=(None, (-102.0, 0.0))):
            phase, expected = self.two_pass(settings)
            transferred = self.process(settings)
            quantitative_phase = self.process(replace(settings, phase_lb_hz=None))
        self.assertEqual(phase.settings.p0_deg, 78.0)
        self.assertEqual(quantitative_phase.settings.p0_deg, -102.0)
        self.assertEqual(transferred.settings.p0_deg, 78.0)
        self.assert_trace_matches(transferred, expected)

    def test_invalid_phase_broadening(self):
        for lb in (nan, inf, -inf, -1):
            with self.subTest(phase_lb_hz=lb):
                with self.assertRaisesRegex(ValueError, "phase_lb_hz must be finite and nonnegative"):
                    self.process(replace(SETTINGS, phase_lb_hz=lb))
        with self.assertRaisesRegex(ValueError, "phase_lb_hz requires an autophase_method"):
            self.process(replace(SETTINGS, phase_lb_hz=0, autophase_method=None))

    def test_linewidth_and_phase_check_on_quantitative_trace(self):
        spectrum = self.process(replace(SETTINGS, phase_lb_hz=1))
        self.assertEqual(spectrum.settings.lb_hz, 0)
        for bounds in REGIONS:
            with self.subTest(bounds=bounds):
                measured = measure_region(spectrum, "fixture peak", bounds, NOISE)
                self.assertIsNone(measured.fwhm_reason)
                self.assertGreater(measured.fwhm_hz, 0)
                phase = check_region_stability(spectrum, bounds, raw=self.raw, noise_bounds_ppm=NOISE).phase
                self.assertIn(phase.status, ("pass", "flag"), phase.reason)
                self.assertIsNotNone(phase.value)
                self.assertTrue(phase.details)

    def test_none_matches_omitted_phase_broadening(self):
        settings = replace(SETTINGS, lb_hz=1)
        default = self.process(settings)
        explicit = self.process(replace(settings, phase_lb_hz=None))
        self.assertIsNone(default.settings.phase_lb_hz)
        np.testing.assert_array_equal(explicit.data, default.data)
        np.testing.assert_array_equal(explicit.smooth_baseline_curve, default.smooth_baseline_curve)
        self.assertEqual(explicit.lineshape_phase_fit, default.lineshape_phase_fit)
        self.assertEqual(explicit.smooth_baseline_fit, default.smooth_baseline_fit)
        # Areas from the default path at 0b459a0, before phase_lb_hz existed.
        np.testing.assert_allclose([integrate(default, bounds).area for bounds in REGIONS],
                                   [0.34859669436599183, 0.25871718794247456], rtol=1e-10, atol=0)
