import json
import unittest
from dataclasses import asdict, replace

import numpy as np

from spinivo import read
from spinivo.filters import resolve_delay
from spinivo.processing import (
    QUANTITATIVE_SETTINGS, ProcessingSettings, _baseline_curve, _smooth_levels, baseline_shape_areas, integrate, process,
)
from tests.integral_fixtures import SETTINGS, make_baseline_raw, make_raw, ppm_bounds
from tests.test_readers import BRUKER, JEOL


# QFY6L9, FID02 Pattern1 0619_2019, matching processed -1-2.jdf (peak 2.46256222 ppm):
# x_offset is 5 ppm in both files; the processed ppm endpoints are centered on 5.
JEOL_REFERENCE_SHIFT_PPM = 0.0
# Zenodo 17083860, 833-1-qnmr/10: pdata/1/procs SF = 400.330037627532 MHz,
# acqus SFO1 = 400.33247204 MHz. Use (SFO1 - SF)*1e6/SF minus raw O1/BF1.
# pdata/1/peaklist.xml reports 1.185527 ppm.
BRUKER_REFERENCE_SHIFT_PPM = -0.09399185873645344


class ProcessingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.jeol = read(JEOL)
        cls.bruker = read(BRUKER)

    def synthetic(self, lines, size=4096, sw=4096.0):
        time = np.arange(size) / sw
        fid = sum(amplitude * np.exp((-4 + 2j * frequency) * np.pi * time) for frequency, amplitude in lines)
        return replace(self.jeol, fid=fid, complex_points=size, acquisition_time_s=size / sw,
                       spectral_width_hz=sw, observe_frequency_mhz=400.0, carrier_ppm=5.0,
                       digital_filter={"digital_filter": False})

    def test_vendor_ppm_axes_and_magnitude_peaks(self):
        for raw, expected, shift, source in (
            (self.jeol, 2.46256222, JEOL_REFERENCE_SHIFT_PPM, "matching processed -1-2.jdf"),
            (self.bruker, 1.185527, BRUKER_REFERENCE_SHIFT_PPM, "pdata/1/procs SF and acqus SFO1"),
        ):
            with self.subTest(vendor=raw.vendor):
                settings = ProcessingSettings(reference_shift_ppm=shift, reference_source=source)
                spectrum = process(raw, settings)
                self.assertTrue(np.all(np.diff(spectrum.ppm) < 0))
                self.assertTrue(np.all(np.diff(spectrum.hz) < 0))
                np.testing.assert_allclose(spectrum.ppm, raw.carrier_ppm + spectrum.hz / raw.observe_frequency_mhz
                                           + shift, rtol=0, atol=5e-15)
                self.assertAlmostEqual(spectrum.ppm[np.argmax(abs(spectrum.data))], expected, delta=0.002)
                self.assertEqual(spectrum.settings, replace(settings, pivot_ppm=raw.carrier_ppm + shift))
                self.assertEqual(spectrum.delay, resolve_delay(raw))
                self.assertEqual(spectrum.acquired_points, raw.complex_points)
                self.assertEqual(spectrum.transformed_points, 2 * raw.complex_points)
                self.assertEqual(len(spectrum.data), spectrum.transformed_points)
                self.assertEqual(spectrum.dwell_s, 1 / raw.spectral_width_hz)
                self.assertEqual(spectrum.first_point_weight, 0.5)
                self.assertEqual(spectrum.intensity_note, raw.intensity_note)

    def test_fixed_area_ratios_linearity_and_zero_fill(self):
        lines = ((-930.3, 1.0), (260.7, 2.0), (1220.2, 3.5))
        raw = self.synthetic(lines)
        bounds = [(5 + (frequency - 75) / 400, 5 + (frequency + 75) / 400) for frequency, _ in lines]
        spectra = (
            process(raw),
            process(replace(raw, fid=2 * raw.fid)),
            process(raw, ProcessingSettings(zero_fill_factor=4)),
        )
        areas, doubled, filled = (np.array([integrate(spectrum, region).area for region in bounds])
                                  for spectrum in spectra)
        np.testing.assert_allclose(areas / areas[0], [amplitude for _, amplitude in lines], rtol=0.005, atol=0)
        np.testing.assert_allclose(doubled, 2 * areas, rtol=1e-12, atol=0)
        np.testing.assert_allclose(filled, areas, rtol=0.001, atol=0)
        self.assertEqual(spectra[2].transformed_points, 4 * raw.complex_points)
        self.assertEqual(spectra[2].settings.zero_fill_factor, 4)

    def test_sampled_band_area_and_apodization(self):
        for size, factor in ((4096, 2), (4095, 1)):
            with self.subTest(size=size, factor=factor):
                raw = self.synthetic(((450.0, 6.0),), size=size)
                original = raw.fid.copy()
                settings = ProcessingSettings(lb_hz=2.5, zero_fill_factor=factor)
                spectrum = process(raw, settings)
                manually_windowed = original * np.exp(-np.pi * settings.lb_hz * np.arange(size) / 4096)
                expected = process(replace(raw, fid=manually_windowed), replace(settings, lb_hz=0))
                np.testing.assert_allclose(spectrum.data, expected.data, rtol=1e-12, atol=1e-12)
                bounds = (spectrum.ppm[-1], spectrum.ppm[0])
                area = integrate(spectrum, bounds).area
                expected_area = 0.5 * original[0].real
                self.assertAlmostEqual(area, expected_area, delta=0.005 * abs(expected_area))
                self.assertAlmostEqual(integrate(expected, bounds).area, expected_area,
                                       delta=0.005 * abs(expected_area))
                self.assertEqual(spectrum.settings.lb_hz, 2.5)
                np.testing.assert_array_equal(raw.fid, original)

    def test_integration_at_lower_grid_edge(self):
        spectrum = process(self.synthetic(((-2047.4, -6.0),)))
        low, inner = spectrum.ppm[-1], spectrum.ppm[-2]
        outside = (spectrum.ppm_limits[0] + low) / 2
        with self.assertRaisesRegex(ValueError, "region reaches the band edge"):
            integrate(spectrum, (outside, inner))

        width = inner - low
        area = integrate(spectrum, (low + width / 4, inner - width / 4)).area
        hz_width = spectrum.hz[-2] - spectrum.hz[-1]
        expected = hz_width / 2 * spectrum.data.real[-2:].mean()
        self.assertLess(expected, 0)
        self.assertAlmostEqual(area, expected)

    def test_manual_phase_preserves_phase_at_referenced_pivot(self):
        raw = self.synthetic(((450.0, 1.0),))
        settings = ProcessingSettings(reference_shift_ppm=0.125, p0_deg=23.0, pivot_ppm=6.25)
        zero_order = process(raw, settings)
        phased = process(raw, replace(settings, p1_deg=77.0))
        at_pivot = np.flatnonzero(phased.ppm == settings.pivot_ppm).item()
        self.assertAlmostEqual(np.angle(phased.data[at_pivot], deg=True), settings.p0_deg, delta=0.001)
        self.assertEqual(phased.data[at_pivot], zero_order.data[at_pivot])
        ramp = np.exp(1j * np.deg2rad(77 * (settings.pivot_ppm - phased.ppm) / (4096 / 400)))
        np.testing.assert_allclose(phased.data, zero_order.data * ramp, atol=1e-12)
        self.assertEqual(phased.settings, replace(settings, p1_deg=77.0))

    def test_baseline_and_signed_area_with_interpolated_bounds(self):
        raw = self.synthetic(((0, 1),), size=1024)
        hz = 2048 - 4 * np.arange(raw.complex_points)
        ppm = 5 + hz / 400
        signal = np.where(abs(hz) < 30, 2.5, 0.0)
        baseline = -7 + 0.5 * ppm
        imaginary = np.sin(hz / 100)
        desired = signal + baseline + 1j * imaginary
        fid = np.fft.fft(np.fft.ifftshift(desired)) * 4
        fid[0] *= 2
        raw = replace(raw, fid=fid)
        settings = ProcessingSettings(zero_fill_factor=1, baseline_windows_ppm=((2, 1), (8, 9)))
        phased = process(raw, settings)
        self.assertEqual(phased.settings.baseline_windows_ppm, settings.baseline_windows_ppm)
        np.testing.assert_allclose(phased.baseline_coefficients, (-7, 0.5), atol=1e-12)
        np.testing.assert_allclose(phased.data, signal + 1j * imaginary, atol=1e-12)

        uncorrected = process(raw, replace(settings, baseline_windows_ppm=()))
        low, high = 1.123, 1.239
        region = integrate(uncorrected, (high, low))
        self.assertEqual(region, integrate(uncorrected, (low, high)))
        self.assertEqual(region.bounds_ppm, (low, high))
        np.testing.assert_allclose(region.bounds_hz, [(p - 5) * 400 for p in (low, high)])
        self.assertEqual(region.unit, "intensity*Hz")
        self.assertAlmostEqual(region.area, (-7 + 0.5 * (low + high) / 2) * (high - low) * 400)

    def test_automatic_baseline_constant_offset(self):
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        bounds = ppm_bounds((-100, 100))
        raw = make_baseline_raw()
        spectrum = process(raw, settings, regions_ppm=(bounds,))
        fit = spectrum.smooth_baseline_fit
        selection = fit.selection
        self.assertEqual(selection.branch, "constant")
        self.assertIsNone(selection.reason)
        self.assertEqual(selection.median_residual, 0.25)
        self.assertLessEqual(selection.residual_mad, 2 * selection.noise)
        self.assertGreaterEqual(selection.adjacent_pairs, 11)
        self.assertEqual(selection.eligible_blocks, len(fit.anchors))
        self.assertEqual(fit.block_size_points, 8)
        self.assertEqual(fit.half_block_selection.branch, "constant")
        np.testing.assert_allclose(spectrum.data.real, np.exp(-0.5 * (spectrum.hz / 2) ** 2), atol=1e-12)

        bridged = process(raw, replace(settings, baseline_bridge=True), regions_ppm=(bounds,))
        self.assertFalse(bridged.smooth_baseline_fit.regions[0].bridged)
        self.assertEqual(bridged.smooth_baseline_fit, fit)
        self.assertEqual(bridged.smooth_baseline_curve, spectrum.smooth_baseline_curve)
        np.testing.assert_array_equal(bridged.data, spectrum.data)

    def test_automatic_baseline_constant_offset_with_lorentzian_wings(self):
        raw = make_raw({"id": "wing-offset"})
        offset = 0.025
        fid = raw.fid.copy()
        # A first-point increment adds a real constant after half weighting and FFT scaling.
        fid[0] += 2 * offset * raw.spectral_width_hz
        bounds = ppm_bounds((-100, 100))
        spectrum = process(replace(raw, fid=fid), SETTINGS, regions_ppm=(bounds,))
        self.assertEqual(len(spectrum.lineshape_phase_fit.peaks), 3)
        selection = spectrum.smooth_baseline_fit.selection
        self.assertEqual(selection.branch, "constant")
        self.assertAlmostEqual(selection.median_residual, offset, delta=3e-6)
        self.assertEqual(selection.eligible_blocks, 708)

    def test_automatic_baseline_curved_background(self):
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        bounds = ppm_bounds((-40, 40))
        spectrum = process(make_baseline_raw(curvature=0.001), settings, regions_ppm=(bounds,))
        selection = spectrum.smooth_baseline_fit.selection
        self.assertEqual(selection.branch, "smooth")
        self.assertIn("Residual MAD", selection.reason)
        self.assertGreaterEqual(selection.adjacent_pairs, 11)
        self.assertGreater(selection.residual_mad, 2 * selection.noise)
        in_window = abs(spectrum.hz) <= 40
        np.testing.assert_allclose(spectrum.data.real[in_window],
                                   np.exp(-0.5 * (spectrum.hz[in_window] / 2) ** 2), atol=1e-5)

    def test_baseline_shape_areas_survive_a_settings_json_round_trip(self):
        settings = ProcessingSettings(**json.loads(json.dumps(asdict(QUANTITATIVE_SETTINGS))))
        self.assertEqual(settings.baseline_windows_ppm, [])
        bounds = (2.90, 3.01)
        original = process(self.jeol, QUANTITATIVE_SETTINGS, regions_ppm=(bounds,))
        rebuilt = process(self.jeol, settings, regions_ppm=(bounds,))
        self.assertIsNotNone(baseline_shape_areas(rebuilt, bounds))
        self.assertEqual(baseline_shape_areas(rebuilt, bounds), baseline_shape_areas(original, bounds))

    def test_baseline_shape_areas_use_the_accepted_anchor_levels(self):
        settings = replace(QUANTITATIVE_SETTINGS, reference_shift_ppm=0.01, reference_source="fixture")
        bounds = (2.91, 3.02)
        spectrum = process(self.jeol, settings, regions_ppm=(bounds,))
        phased = process(self.jeol, replace(settings, smooth_baseline=False), regions_ppm=(bounds,))
        self.assertEqual(spectrum.settings.p0_deg, phased.settings.p0_deg)
        self.assertEqual(spectrum.settings.p1_deg, phased.settings.p1_deg)
        ppm, levels = np.asarray(spectrum.smooth_baseline_fit.anchors).T
        lower, upper = np.where(ppm < bounds[0])[0][0], np.where(ppm > bounds[1])[0][-1]
        line = np.interp(spectrum.ppm, ppm[[lower, upper]], levels[[lower, upper]])
        areas = baseline_shape_areas(spectrum, bounds[::-1])
        self.assertEqual(areas.production, integrate(spectrum, bounds).area)
        for area, curve in ((areas.anchor_median, np.median(levels)), (areas.bracketing_line, line)):
            expected = integrate(replace(phased, data=phased.data - curve), bounds).area
            np.testing.assert_allclose(area, expected, rtol=1e-12)
        self.assertIsNone(baseline_shape_areas(phased, bounds))

        edge = (spectrum.ppm[-1], spectrum.ppm[-1] + 0.02)
        edge_spectrum = process(self.jeol, settings, regions_ppm=(bounds, edge))
        edge_areas = baseline_shape_areas(edge_spectrum, edge)
        self.assertIsNotNone(edge_areas.anchor_median)
        self.assertIsNone(edge_areas.bracketing_line)

    def test_automatic_baseline_bridge_uses_close_anchors(self):
        raw = make_baseline_raw(curvature=0.001)
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True, baseline_bridge=False)
        bounds = ppm_bounds((-40, 40))
        original = process(raw, settings, regions_ppm=(bounds,))
        fit = original.smooth_baseline_fit
        self.assertFalse(fit.regions[0].bridged)

        # Reconstruct the legacy fit at block centers before sampling it on the ppm grid.
        blocks = original.ppm.reshape(-1, fit.block_size_points)
        centers = (blocks[:, 0] + blocks[:, -1]) / 2
        anchor_ppm, anchor_levels = np.asarray(fit.anchors).T
        anchors = np.isin(centers, anchor_ppm)
        levels = np.zeros_like(centers)
        levels[anchors] = anchor_levels
        smooth = _smooth_levels(levels, anchors, (0.025 / (centers[0] - centers[1])) ** 4)
        legacy = np.interp(original.ppm, centers[::-1], smooth[::-1])
        np.testing.assert_array_equal(original.smooth_baseline_curve.view(np.uint64), legacy.view(np.uint64))
        uncorrected = process(raw, replace(settings, smooth_baseline=False), regions_ppm=(bounds,))
        np.testing.assert_array_equal(original.data.view(np.uint64), (uncorrected.data - legacy).view(np.uint64))

        spectrum = process(raw, replace(settings, baseline_bridge=True), regions_ppm=(bounds,))
        fit = spectrum.smooth_baseline_fit
        region, = fit.regions
        self.assertEqual(fit.selection.branch, "smooth")
        self.assertTrue(region.bridged)
        low, high = bounds
        lower = np.flatnonzero(anchors & (centers < low))[0]
        upper = np.flatnonzero(anchors & (centers > high))[-1]
        between = (centers[lower] <= spectrum.ppm) & (spectrum.ppm <= centers[upper])
        expected = np.interp(spectrum.ppm[between], centers[[lower, upper]], smooth[[lower, upper]])
        self.assertGreater(np.max(abs(legacy[between] - expected)), 1e-6)
        np.testing.assert_allclose(spectrum.smooth_baseline_curve[between], expected, atol=1e-15, rtol=0)
        np.testing.assert_array_equal(spectrum.smooth_baseline_curve[~between], legacy[~between])
        self.assertAlmostEqual(region.half_block_area, 5.013121470831763, delta=1e-9)
        self.assertAlmostEqual(region.recovery.value, 0.997772781239347, delta=1e-9)
        self.assertEqual(region.recovery.placement_sensitivity, 0.0)
        nominal = integrate(original, bounds).area
        straight = integrate(spectrum, bounds).area
        self.assertEqual(original.smooth_baseline_fit.regions[0].recovery.placement_sensitivity,
                         abs(straight - nominal) / nominal)

    def test_automatic_baseline_bridge_skips_distant_anchor(self):
        raw = make_baseline_raw(curvature=0.001)
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        bounds = ppm_bounds((-1, 20))
        original = process(raw, settings, regions_ppm=(bounds,))
        spectrum = process(raw, replace(settings, baseline_bridge=True), regions_ppm=(bounds,))
        fit = spectrum.smooth_baseline_fit
        self.assertEqual(fit.selection.branch, "smooth")
        region, = fit.regions
        below, above = region.anchor_distances_ppm
        half_width = (bounds[1] - bounds[0]) / 2
        self.assertLessEqual(above, half_width)
        self.assertGreater(below, half_width)
        self.assertFalse(region.bridged)
        np.testing.assert_array_equal(spectrum.smooth_baseline_curve, original.smooth_baseline_curve)
        np.testing.assert_array_equal(spectrum.data, original.data)
        self.assertAlmostEqual(region.recovery.placement_sensitivity, 2.9208705163430437e-6, delta=1e-10)

    def test_baseline_bridge_settings(self):
        self.assertTrue(QUANTITATIVE_SETTINGS.baseline_bridge)
        self.assertFalse(ProcessingSettings().baseline_bridge)
        self.assertFalse(SETTINGS.baseline_bridge)

    def test_spike_anchor_lock_settings(self):
        self.assertTrue(QUANTITATIVE_SETTINGS.lock_spike_anchors)
        self.assertFalse(ProcessingSettings().lock_spike_anchors)
        self.assertFalse(SETTINGS.lock_spike_anchors)

    def test_locked_spike_refit_reports_the_nominal_branch_it_uses(self):
        spectrum = process(make_baseline_raw(), ProcessingSettings(zero_fill_factor=1))
        regions = (ppm_bounds((-100, 100)),)
        nominal = _baseline_curve(spectrum, regions, 8)
        self.assertEqual(nominal.selection.branch, "constant")
        trial = replace(spectrum, data=spectrum.data + 0.1 / (1 + spectrum.hz ** 2))
        free = _baseline_curve(trial, regions, 8)
        self.assertEqual(free.selection.branch, "smooth")
        locked = _baseline_curve(trial, regions, 8, anchor_lock=(nominal.anchor_mask, nominal.selection.branch))
        self.assertEqual(locked.selection.branch, "constant")
        self.assertGreater(locked.curve, nominal.curve)

    def test_automatic_baseline_with_too_few_adjacent_anchors(self):
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        bounds = ppm_bounds((-2030, 2030))
        spectrum = process(make_baseline_raw(), settings, regions_ppm=(bounds,))
        selection = spectrum.smooth_baseline_fit.selection
        self.assertEqual(selection.branch, "smooth")
        self.assertIn("Fewer than 11 adjacent eligible block pairs", selection.reason)
        self.assertEqual((selection.eligible_blocks, selection.adjacent_pairs), (4, 2))
        self.assertAlmostEqual(integrate(spectrum, bounds).area, np.sqrt(2 * np.pi) * 2, places=5)

    def test_lineshape_autophase_recovers_phase_with_curved_baseline(self):
        lines = ((-1790.3, 1.0), (555.8, 2.0), (1030.2, 3.0), (1740.6, 1.5))
        raw = self.synthetic(lines)
        settings = ProcessingSettings(zero_fill_factor=1, autophase_method="lineshape", pivot_ppm=4.25)
        clean = process(raw, replace(settings, autophase_method=None))
        p0, p1 = -146.0, 334.0
        position = (settings.pivot_ppm - clean.ppm) / (4096 / 400)
        desired = clean.data * np.exp(-1j * np.deg2rad(p0 + p1 * position))
        desired += (0.015 + 0.007j) * (clean.hz / 1024) ** 2
        fid = np.fft.fft(np.fft.ifftshift(desired)) / (clean.transformed_points * clean.dwell_s)
        fid[0] *= 2

        phased = process(replace(raw, fid=fid), settings)
        self.assertAlmostEqual(phased.settings.p0_deg, p0, delta=1)
        self.assertAlmostEqual(phased.settings.p1_deg, p1, delta=1)
        fit = phased.lineshape_phase_fit
        self.assertEqual(len(fit.peaks), 4)
        np.testing.assert_allclose([peak.ppm for peak in fit.peaks],
                                   [5 + hz / 400 for hz, _ in reversed(lines)], rtol=0, atol=0.001)
        self.assertEqual(fit.slope_deg, phased.settings.p1_deg)
        self.assertLess(fit.phase_rms_deg, 1)

    def test_lineshape_rejects_ambiguous_distant_peaks(self):
        raw = self.synthetic(((-1700, np.exp(2j)), (1700, np.exp(-2j))))
        with self.assertRaisesRegex(ValueError, "ambiguous phase"):
            process(raw, ProcessingSettings(autophase_method="lineshape"), regions_ppm=((4.9, 5.1),))

    def test_lineshape_single_peak_requires_nearby_region_edges(self):
        raw = self.synthetic(((450.0, np.exp(0.4j)),))
        settings = ProcessingSettings(autophase_method="lineshape")
        spectrum = process(raw, settings, regions_ppm=((6.10, 6.15),))
        fit = spectrum.lineshape_phase_fit
        self.assertEqual(len(fit.peaks), 1)
        self.assertEqual(fit.method, "single_peak")
        self.assertEqual(fit.slope_deg, 0.0)
        self.assertIsNone(fit.phase_rms_deg)
        self.assertEqual(fit.note, "slope not determined")
        self.assertEqual(spectrum.settings.p1_deg, 0)
        self.assertAlmostEqual(spectrum.settings.p0_deg, -np.rad2deg(0.4), delta=1)
        with self.assertRaisesRegex(ValueError, "insufficient peak span"):
            process(raw, settings, regions_ppm=((6.10, 6.30),))

    def test_lineshape_fits_isolated_peaks_beyond_first_32_candidates(self):
        doublets = [(center + offset, 4.0) for center in np.linspace(-1300, 1300, 34) for offset in (-2, 2)]
        isolated = [(-1700, 1.0), (1700, 1.0)]
        spectrum = process(self.synthetic(doublets + isolated), ProcessingSettings(autophase_method="lineshape"))
        np.testing.assert_allclose([peak.ppm for peak in spectrum.lineshape_phase_fit.peaks],
                                   [9.25, 0.75], rtol=0, atol=0.001)

    def test_lineshape_on_spectrum_with_fewer_than_100_points(self):
        raw = self.synthetic(((-25.4, np.exp(0.4j)), (24.2, 1.2 * np.exp(0.4j))), size=48, sw=96)
        spectrum = process(raw, ProcessingSettings(autophase_method="lineshape"))
        self.assertEqual(spectrum.transformed_points, 96)
        self.assertEqual(len(spectrum.lineshape_phase_fit.peaks), 2)
        self.assertAlmostEqual(spectrum.settings.p0_deg, -np.rad2deg(0.4), delta=1)
        self.assertAlmostEqual(spectrum.settings.p1_deg, 0, delta=1)

    def test_autophase_records_applied_phase_about_pivot(self):
        raw = self.synthetic(((-1536, np.exp(0.4j)), (768, 2 * np.exp(0.5j))), size=1024)
        settings = ProcessingSettings(autophase_method="acme", pivot_ppm=4.0)
        phased = process(raw, settings)
        unphased = process(raw)
        applied = phased.settings
        self.assertEqual(applied.autophase_method, settings.autophase_method)
        self.assertEqual(applied.pivot_ppm, settings.pivot_ppm)
        self.assertGreater(abs(applied.p0_deg), 1)
        ramp = applied.p0_deg + applied.p1_deg * (applied.pivot_ppm - phased.ppm) / (4096 / 400)
        np.testing.assert_allclose(phased.data, unphased.data * np.exp(1j * np.deg2rad(ramp)), atol=1e-12)

    def test_autophase_preserves_area_scaling(self):
        lines = ((-1536, 0.4 * np.exp(0.4j)), (768, 0.8 * np.exp(0.5j)))
        raw = self.synthetic(lines, size=1024)
        settings = ProcessingSettings(autophase_method="acme")
        original = process(raw, settings)
        doubled = process(replace(raw, fid=2 * raw.fid), settings)
        for frequency, _ in lines:
            with self.subTest(frequency=frequency):
                bounds = (5 + (frequency - 75) / 400, 5 + (frequency + 75) / 400)
                expected = 2 * integrate(original, bounds).area
                self.assertGreater(expected, 0)
                self.assertAlmostEqual(integrate(doubled, bounds).area, expected, delta=0.005 * expected)

    def test_autophase_corrects_inverted_line_and_records_phase(self):
        raw = self.synthetic(((450, 1),))
        inverted = replace(raw, fid=-raw.fid)
        settings = ProcessingSettings(autophase_method="acme")
        upright = process(raw, settings)
        corrected = process(inverted, settings)
        bounds = (5 + (450 - 75) / 400, 5 + (450 + 75) / 400)
        area = integrate(corrected, bounds).area
        expected = integrate(upright, bounds).area
        self.assertGreater(area, 0)
        self.assertAlmostEqual(area, expected, delta=0.005 * expected)
        applied = replace(corrected.settings, autophase_method=None)
        np.testing.assert_allclose(corrected.data, process(inverted, applied).data, atol=1e-12)

    def test_refuses_unresolved_delay(self):
        with self.assertRaisesRegex(ValueError, "unresolved digital filter delay"):
            process(replace(self.bruker, digital_filter=None))

    def test_refuses_unusable_baseline_limit(self):
        raw = self.synthetic(((0, 1),))
        for limit in (float("nan"), float("inf"), -0.01, None, True):
            with self.subTest(limit=limit), self.assertRaisesRegex(ValueError, "baseline_limit must be finite"):
                process(raw, baseline_limit=limit)
        self.assertEqual(process(raw, baseline_limit=0).baseline_limit, 0)

    def test_refuses_missing_reference_metadata_and_bounds_outside_band(self):
        raw = self.synthetic(((0, 1),))
        with self.assertRaisesRegex(ValueError, "carrier_ppm"):
            process(replace(raw, carrier_ppm=None))
        spectrum = process(raw)
        with self.assertRaisesRegex(ValueError, "region reaches the band edge"):
            integrate(spectrum, (spectrum.ppm_limits[0] - 0.1, spectrum.ppm_limits[1]))
