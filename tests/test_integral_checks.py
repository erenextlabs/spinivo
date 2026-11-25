import unittest
from dataclasses import replace
from math import inf, nan, nextafter
from unittest.mock import patch

import numpy as np

from spinivo.integral_checks import (
    PILOT, LockedBaselineRecovery, _measure_baseline_recovery, check_baseline_recovery, check_region_stability,
)
from spinivo.processing import ProcessingSettings, integrate, process
from tests.integral_fixtures import NOISE_BOUNDS_HZ, SETTINGS, make_baseline_raw, make_raw, ppm_bounds


class BaselineRecoveryTests(unittest.TestCase):
    def test_clean_singlet_recovers_added_area(self):
        bounds = ppm_bounds((-100, 100))
        spectrum = process(make_raw({"id": "clean"}), SETTINGS, regions_ppm=(bounds,))
        recovery = check_baseline_recovery(spectrum, bounds[::-1])
        self.assertEqual(recovery, spectrum.smooth_baseline_fit.regions[0].recovery)
        self.assertAlmostEqual(recovery.value, 1, delta=0.001)
        self.assertAlmostEqual(recovery.recovery_25, 1, delta=0.002)
        self.assertAlmostEqual(recovery.recovery_75, 1, delta=0.002)
        self.assertEqual(recovery.placement_sensitivity, 0.0)
        self.assertEqual(recovery.limit, 0.017)
        self.assertIsNone(recovery.flag)
        self.assertIsNone(recovery.reason)

    def test_shared_recovery_loss_does_not_flag(self):
        bounds = ppm_bounds((-20, 20))
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        spectrum = process(make_baseline_raw(curvature=0.25), settings, regions_ppm=(bounds,))
        recovery = check_baseline_recovery(spectrum, bounds)
        self.assertLess(max(recovery.value, recovery.recovery_25, recovery.recovery_75), 0.99)
        self.assertLess(recovery.position_sensitivity, recovery.limit)
        self.assertLess(recovery.placement_sensitivity, recovery.limit)
        self.assertIsNone(recovery.flag)
        self.assertIsNone(recovery.reason)

    def test_position_sensitivity_flags_with_zero_placement_change(self):
        bounds = ppm_bounds((-15, 25))
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True, baseline_bridge=True)
        spectrum = process(make_baseline_raw(curvature=0.25), settings, regions_ppm=(bounds,))
        recovery = check_baseline_recovery(spectrum, bounds)
        self.assertEqual(recovery.position_sensitivity, abs(recovery.recovery_25 - recovery.recovery_75))
        self.assertGreater(recovery.position_sensitivity, recovery.limit)
        self.assertEqual(recovery.placement_sensitivity, 0.0)
        self.assertEqual(recovery.flag, "baseline_sensitive")
        self.assertIsNone(recovery.reason)

    def test_spike_tail_still_flags_with_locked_anchors(self):
        bounds = ppm_bounds((-30, 7))
        raw = make_baseline_raw(curvature=0.25)
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True, baseline_bridge=True)
        for locked in (False, True):
            with self.subTest(locked=locked):
                spectrum = process(raw, replace(settings, lock_spike_anchors=locked), regions_ppm=(bounds,))
                baseline = check_region_stability(spectrum, bounds).baseline
                recovery = baseline.details.recovery
                self.assertEqual(isinstance(recovery, LockedBaselineRecovery), locked)
                self.assertGreater(recovery.position_sensitivity, 0.017)
                self.assertEqual(recovery.placement_sensitivity, 0.0)
                self.assertLess(baseline.value, baseline.limit)
                self.assertEqual(recovery.flag, "baseline_sensitive")
                self.assertEqual(baseline.status, "flag")

    def test_placement_sensitivity_flags_with_small_position_change(self):
        bounds = ppm_bounds((-200, 200))
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        spectrum = process(make_baseline_raw(curvature=0.02), settings, regions_ppm=(bounds,))
        recovery = check_baseline_recovery(spectrum, bounds)
        self.assertEqual(spectrum.smooth_baseline_fit.selection.branch, "smooth")
        self.assertLess(recovery.position_sensitivity, recovery.limit)
        self.assertGreater(recovery.placement_sensitivity, recovery.limit)
        self.assertEqual(recovery.flag, "baseline_sensitive")
        self.assertIsNone(recovery.reason)

    def test_sensitivity_limit_is_strict(self):
        bounds = ppm_bounds((-100, 100))
        for position, flag in ((0.017, None), (nextafter(0.017, inf), "baseline_sensitive")):
            with self.subTest(position=position), patch("spinivo.integral_checks._recovery_trial",
                                                       side_effect=(0.98, 0.0, position)):
                recovery = _measure_baseline_recovery(None, bounds, 1.0, None, 1.0)
                self.assertEqual(recovery.position_sensitivity, position)
                self.assertEqual(recovery.limit, 0.017)
                self.assertEqual(recovery.flag, flag)
        with patch("spinivo.integral_checks._recovery_trial", side_effect=(0.98, 0.0, 0.0)):
            recovery = _measure_baseline_recovery(None, bounds, 1000.0, None, 983.0)
        self.assertEqual(recovery.placement_sensitivity, 0.017)
        self.assertIsNone(recovery.flag)

    def test_failed_trial_leaves_a_gap_and_retains_placement_flag(self):
        bounds = ppm_bounds((-100, 100))
        with patch("spinivo.integral_checks._recovery_trial",
                   side_effect=(0.98, ValueError("No quiet anchor blocks"), 0.99)):
            recovery = _measure_baseline_recovery(None, bounds, 1.0, None, 1.03)
        self.assertEqual(recovery.value, 0.98)
        self.assertIsNone(recovery.recovery_25)
        self.assertEqual(recovery.recovery_75, 0.99)
        self.assertIsNone(recovery.position_sensitivity)
        self.assertIn("25% recovery: No quiet anchor blocks", recovery.reason)
        self.assertEqual(recovery.flag, "baseline_sensitive")

    def test_missing_anchor_bracket_leaves_placement_gap(self):
        bounds = ppm_bounds((-2047, 100))
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        spectrum = process(make_baseline_raw(curvature=0.001), settings, regions_ppm=(bounds,))
        recovery = check_baseline_recovery(spectrum, bounds)
        self.assertLess(recovery.position_sensitivity, recovery.limit)
        self.assertIsNone(recovery.placement_sensitivity)
        self.assertIsNone(recovery.flag)
        self.assertIn("Missing baseline anchor bracket", recovery.reason)
        self.assertEqual(check_region_stability(spectrum, bounds).baseline.status, "cannot be assessed")

    def test_explicit_baseline_windows_repeat_the_linear_fit(self):
        bounds = ppm_bounds((-100, 100))
        windows = (ppm_bounds((-300, -200)), ppm_bounds((200, 300)))
        settings = ProcessingSettings(zero_fill_factor=1, baseline_windows_ppm=windows)
        spectrum = process(make_baseline_raw(), settings, regions_ppm=(bounds,))
        recovery = check_baseline_recovery(spectrum, bounds)
        self.assertAlmostEqual(recovery.value, 1, delta=0.002)
        self.assertLess(recovery.value, 1)
        self.assertEqual(recovery.placement_sensitivity, 0.0)
        self.assertIsNotNone(recovery.position_sensitivity)
        self.assertIsNone(recovery.flag)
        self.assertIsNone(recovery.reason)

    def test_unavailable_recovery_reports_reason(self):
        bounds = ppm_bounds((-100, 100))
        raw = make_baseline_raw()
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        cases = (
            (process(raw), "no baseline correction"),
            (process(raw, settings), "No recovery recorded"),
            (process(replace(raw, fid=-raw.fid), settings, regions_ppm=(bounds,)), "finite and positive"),
        )
        for spectrum, reason in cases:
            with self.subTest(reason=reason):
                recovery = check_baseline_recovery(spectrum, bounds)
                self.assertIsNone(recovery.value)
                self.assertIsNone(recovery.flag)
                self.assertEqual(recovery.limit, 0.017)
                self.assertIn(reason, recovery.reason)


class RegionStabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bounds = ppm_bounds((-100, 100))
        cls.noise = ppm_bounds(NOISE_BOUNDS_HZ)
        cls.raw = make_raw({"id": "clean"})
        cls.spectrum = process(cls.raw, SETTINGS, regions_ppm=(cls.bounds,))

    def check(self, spectrum=None, **options):
        options = dict(raw=self.raw, noise_bounds_ppm=self.noise) | options
        return check_region_stability(spectrum if spectrum is not None else self.spectrum, self.bounds, **options)

    def assert_outcome(self, observation, status, reason):
        self.assertEqual(observation.status, status)
        self.assertIn(reason, observation.reason)

    def test_clean_region_passes_without_changing_the_trace(self):
        trace = self.spectrum.data.copy()
        expected = integrate(self.spectrum, self.bounds).area
        result = self.check()
        for name in ("area", "baseline", "bounds", "phase", "snr"):
            item = getattr(result, name)
            self.assertEqual((item.check, item.status, item.threshold_version), (name, "pass", "pilot-2"))
        self.assertEqual((result.nominal_area, result.area.value), (expected, expected))
        alternate = self.spectrum.smooth_baseline_fit.regions[0].half_block_area
        self.assertEqual(result.baseline.value, abs(alternate - expected) / abs(expected))
        np.testing.assert_array_equal(self.spectrum.data, trace)

    def test_metric_limits(self):
        with self.assertRaisesRegex(ValueError, "max_relative_change"):
            replace(PILOT, max_relative_change=nan)
        nominal = self.check()
        for item in (nominal.baseline, nominal.bounds, nominal.phase, nominal.snr):
            setting, direction = ("min_snr", inf) if item.check == "snr" else ("max_relative_change", -inf)
            for limit, status in ((item.value, "pass"), (nextafter(item.value, direction), "flag")):
                thresholds = replace(PILOT, version="test", **{setting: limit})
                result = getattr(self.check(thresholds=thresholds), item.check)
                self.assertEqual((result.status, result.limit, result.threshold_version), (status, limit, "test"))

    def test_four_independent_edge_moves(self):
        result = self.check()
        expected = {"low_out": (-110, 100), "low_in": (-90, 100), "high_in": (-100, 90), "high_out": (-100, 110)}
        self.assertEqual(tuple(trial.move for trial in result.bounds.details), tuple(expected))
        for trial in result.bounds.details:
            np.testing.assert_allclose(trial.bounds_ppm, ppm_bounds(expected[trial.move]), rtol=0, atol=1e-14)
            area = integrate(self.spectrum, trial.bounds_ppm).area
            self.assertEqual(trial.relative_change, (area - result.nominal_area) / abs(result.nominal_area))
        self.assertEqual(result.bounds.value, max(abs(t.relative_change) for t in result.bounds.details))

    def test_incomplete_edges_retain_admissible_flags(self):
        spectrum = replace(self.spectrum, data=np.ones_like(self.spectrum.data))
        high = spectrum.ppm[0]
        for limit, status in ((0.01, "flag"), (0.10, "cannot be assessed")):
            result = check_region_stability(spectrum, (high - 1, high),
                                            thresholds=replace(PILOT, max_relative_change=limit)).bounds
            self.assert_outcome(result, status, "outside")
            self.assertIsNone(result.details[-1].relative_change)
            self.assertIn("outside", result.details[-1].reason)
        result = check_region_stability(spectrum, ppm_bounds((0, 0.0625))).bounds
        self.assertTrue(all("collapse" in trial.reason for trial in result.details[1:3]))
        self.assertAlmostEqual(result.details[0].bounds_ppm[0], ppm_bounds((-0.125, 0))[0])

    def test_missing_baseline_anchors_and_alternate(self):
        fit = self.spectrum.smooth_baseline_fit
        region = fit.regions[0]
        left, right = region.anchor_distances_ppm
        for distances, alternate, reason in (((None, right), region.half_block_area, "left baseline anchor"),
                                             ((left, None), region.half_block_area, "right baseline anchor"),
                                             ((left, right), None, "half-block area")):
            for status, recovery in (("cannot be assessed", region.recovery),
                                     ("flag", replace(region.recovery, position_sensitivity=0.02,
                                                      flag="baseline_sensitive"))):
                changed = replace(region, anchor_distances_ppm=distances, half_block_area=alternate, recovery=recovery)
                spectrum = replace(self.spectrum, smooth_baseline_fit=replace(fit, regions=(changed,)))
                result = self.check(spectrum).baseline
                self.assert_outcome(result, status, reason)
                self.assertEqual((result.details.anchor_distances_ppm, result.details.recovery), (distances, recovery))

    def test_invalid_area_disables_relative_metrics(self):
        for value, reason in ((0, "nonpositive"), (-1, "nonpositive"), (nan, "nonfinite")):
            spectrum = replace(self.spectrum, data=np.full_like(self.spectrum.data, value))
            result = self.check(spectrum)
            np.testing.assert_equal(result.nominal_area, integrate(spectrum, self.bounds).area)
            for item in (result.area, result.baseline, result.bounds, result.phase):
                self.assert_outcome(item, "cannot be assessed", reason)

    def test_injected_phase_error_is_corrected_in_the_right_direction(self):
        raw = make_raw({"id": "phase", "phase_deg": 20})
        windows = (ppm_bounds((-300, -200)), ppm_bounds((200, 300)))
        for requested in (SETTINGS, replace(SETTINGS, smooth_baseline=False, baseline_windows_ppm=windows)):
            spectrum = process(raw, requested, regions_ppm=(self.bounds,))
            result = self.check(spectrum)
            trial, = result.phase.details
            settings = spectrum.settings
            applied = settings.p0_deg + settings.p1_deg * (settings.pivot_ppm - trial.peak_ppm) / (
                spectrum.spectral_width_hz / spectrum.observe_frequency_mhz)
            self.assertAlmostEqual(trial.residual_deg, -20 - applied, delta=0.1)
            self.assertEqual(result.phase.status, "flag")
            self.assertGreater(trial.relative_change, 0.01)
            baseline = (spectrum.smooth_baseline_curve if requested.smooth_baseline
                        else np.polynomial.polynomial.polyval(spectrum.ppm, spectrum.baseline_coefficients))
            restored = (spectrum.data + baseline) * np.exp(1j * np.deg2rad(trial.residual_deg)) - baseline
            expected = integrate(replace(spectrum, data=restored), self.bounds).area
            self.assertAlmostEqual(result.nominal_area * (1 + trial.relative_change), expected, places=13)

    def test_phase_requires_local_fit_evidence(self):
        fit = self.spectrum.lineshape_phase_fit
        for alternative, reason in ((None, "No lineshape"),
                                    (replace(fit, method="single_peak"), "single_peak"),
                                    (replace(fit, peaks=()), "No accepted phase peak")):
            result = self.check(replace(self.spectrum, lineshape_phase_fit=alternative)).phase
            self.assert_outcome(result, "cannot be assessed", reason)

    def test_snr_preserves_acquisition_status_and_reason(self):
        raw = make_raw({"id": "weak", "snr": 100})
        spectrum = process(raw, SETTINGS, regions_ppm=(self.bounds,))
        result = self.check(spectrum, raw=raw).snr
        self.assertLess(result.value, 150)
        self.assertEqual((result.status, result.details.status, result.reason), ("flag", "flag", result.details.reason))

    def test_snr_requires_raw_source_when_spectrum_source_is_known(self):
        raw = replace(self.raw, source=None)
        result = self.check(raw=raw).snr
        self.assert_outcome(result, "cannot be assessed", "Raw source is missing")
        self.assertIsNone(result.value)
        self.assertEqual(self.check(replace(self.spectrum, source=None), raw=raw).snr.status, "pass")

    def test_snr_requires_raw_and_qualified_noise_window(self):
        foreign = replace(self.raw, complex_points=2 * self.raw.complex_points)
        for options, reason in ((dict(raw=None), "raw acquisition"), (dict(raw=foreign), "complex_points"),
                                (dict(noise_bounds_ppm=None), "noise window"),
                                (dict(noise_bounds_ppm=ppm_bounds((300, 302))), "32 acquired points")):
            self.assert_outcome(self.check(**options).snr, "cannot be assessed", reason)
        # 300 to 304 Hz spans exactly 32 acquired points; ppm conversion leaves roundoff below 32.
        self.assertNotEqual(self.check(noise_bounds_ppm=ppm_bounds((300, 304))).snr.status, "cannot be assessed")
        self.assert_outcome(self.check(noise_bounds_ppm=ppm_bounds((300, 303.99))).snr, "cannot be assessed",
                            "32 acquired points")
        # At 512 MHz each acquired point spans 1/4096 ppm.
        raw = replace(self.raw, observe_frequency_mhz=512.0)
        spectrum = process(raw, SETTINGS, regions_ppm=(self.bounds,))
        for points, minimum, status in ((32, 32, "pass"), (32, 33, "cannot be assessed"), (33 - 1e-8, 33, "pass")):
            result = self.check(spectrum, raw=raw, thresholds=replace(PILOT, min_noise_points=minimum),
                                noise_bounds_ppm=(6, 6 + points / 4096)).snr
            self.assertEqual(result.status, status)
