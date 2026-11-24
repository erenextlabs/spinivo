import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from spinivo import RawFID
from spinivo.acquisition import SignalEvidence, SignalMeasurement, check_resolution, check_snr
from spinivo.measurement import measure_region
from spinivo.processing import ProcessingSettings, integrate, process


class MeasurementTests(unittest.TestCase):
    def setUp(self):
        self.rate = 24.0
        self.amplitude = 3.0
        size = 4096
        sw = 4096.0
        frequency = 200.0
        observe = 400.0
        carrier = 5.0
        time = np.arange(size) / sw
        fid = self.amplitude * np.exp((-self.rate + 2j * np.pi * frequency) * time)
        self.raw = RawFID(
            vendor="jeol", source=Path("synthetic"), nucleus="1H", observe_frequency_mhz=observe,
            spectral_width_hz=sw, carrier_ppm=carrier, complex_points=size, acquisition_time_s=size / sw,
            relaxation_delay_s=None, pulse_width_us=None, flip_angle_deg=None, pulse_90_us=None,
            scans=None, dummy_scans=None, receiver_gain=None, temperature_k=None, pulse_program=None,
            digital_filter={"digital_filter": False}, intensity_note="synthetic", fid=fid,
            source_fields={}, missing=(),
        )
        self.center = carrier + frequency / observe
        self.signal_bounds = (5.4, 5.6)
        self.noise_bounds = (0.5, 2.5)

    def test_lorentzian_width_position_and_signed_area(self):
        spectrum = process(self.raw)
        result = measure_region(spectrum, "analyte", self.signal_bounds[::-1], self.noise_bounds[::-1])
        integral = integrate(spectrum, self.signal_bounds)
        width = self.rate / np.pi
        self.assertEqual(result.name, "analyte")
        self.assertEqual(result.signal_bounds_ppm, self.signal_bounds)
        self.assertEqual(result.noise_bounds_ppm, self.noise_bounds)
        self.assertEqual(result.peak_ppm, self.center)
        self.assertEqual(result.area, integral.area)
        self.assertEqual(result.area_unit, integral.unit)
        self.assertAlmostEqual(result.fwhm_hz, width, delta=0.02 * width)
        self.assertIsNone(result.fwhm_reason)

    def test_smooth_baseline_removes_curvature_and_preserves_area_ratio(self):
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True)
        grid = process(self.raw, replace(settings, smooth_baseline=False))
        peaks = 3 * np.exp(-0.5 * ((grid.hz - 200) / 4) ** 2) + 7.5 * np.exp(-0.5 * ((grid.hz - 1200) / 4) ** 2)
        desired = peaks + 0.25 + 0.002 * (grid.ppm - 5) ** 2
        fid = np.fft.fft(np.fft.ifftshift(desired)) / (grid.transformed_points * grid.dwell_s)
        fid[0] *= 2
        regions = (self.signal_bounds, (7.9, 8.1))
        spectrum = process(replace(self.raw, fid=fid), settings, regions_ppm=regions)
        first, second = [measure_region(spectrum, "line", bounds, self.noise_bounds) for bounds in regions]
        self.assertAlmostEqual(second.area / first.area, 2.5, delta=0.005 * 2.5)
        self.assertAlmostEqual(first.height, 3.0, delta=0.001)
        np.testing.assert_allclose(spectrum.data.real, peaks, atol=0.001, rtol=0)
        fit = spectrum.smooth_baseline_fit
        self.assertFalse(any(low <= center <= high for center, _ in fit.anchors for low, high in regions))
        for result, diagnostic in zip((first, second), fit.regions):
            self.assertEqual(result.area, integrate(spectrum, result.signal_bounds_ppm).area)
            self.assertEqual(result.baseline_region, diagnostic)
            self.assertTrue(all(distance > 0 for distance in diagnostic.anchor_distances_ppm))
            self.assertLess(abs(diagnostic.half_block_relative_change), 0.005)

    def test_height_noise_and_acquisition_checks(self):
        size = self.raw.complex_points
        rng = np.random.default_rng(42)
        sigma = 0.01
        fid = self.raw.fid + sigma * (rng.standard_normal(size) + 1j * rng.standard_normal(size))
        raw = replace(self.raw, fid=fid)
        spectrum = process(raw)
        result = measure_region(spectrum, "analyte", self.signal_bounds, self.noise_bounds)
        height = self.amplitude / self.rate
        # Half weighting the first point gives noise variance proportional to N - 3/4.
        noise_sd = sigma * np.sqrt(size - 0.75) * spectrum.dwell_s
        self.assertAlmostEqual(result.height, height, delta=0.1 * height)
        self.assertAlmostEqual(result.noise_sd, noise_sd, delta=0.1 * noise_sd)
        low, high = self.noise_bounds
        noise = spectrum.data.real[(low <= spectrum.ppm) & (spectrum.ppm <= high)]
        self.assertAlmostEqual(result.noise_sd, noise.std(ddof=1), places=16)

        measurement = result.to_signal_measurement()
        self.assertEqual(measurement, SignalMeasurement(
            height=result.height, noise_sd=result.noise_sd, fwhm_hz=result.fwhm_hz,
            source=result.source, noise_bounds_ppm=self.noise_bounds,
            transformed_points=spectrum.transformed_points,
        ))
        signal = SignalEvidence(result.name, measurement=measurement)
        self.assertEqual(check_snr(raw, signal).status, "pass")
        self.assertEqual(check_resolution(raw, signal).status, "pass")

    def test_unavailable_noise_keeps_signal_measurements(self):
        spectrum = process(self.raw)
        for bounds, reason in (((5.5, 5.7), "overlaps"), ((5.6, 5.7), "overlaps"),
                               ((1.9999, 2.0001), "fewer than two"), ((-0.5, 2.5), "outside")):
            with self.subTest(bounds=bounds):
                result = measure_region(spectrum, "analyte", self.signal_bounds, bounds)
                self.assertIsNone(result.noise_sd)
                self.assertIn(reason, result.noise_sd_reason)
                self.assertGreater(result.height, 0)
                self.assertGreater(result.fwhm_hz, 0)

    def test_absent_noise_window_keeps_area_linewidth_and_provenance(self):
        spectrum = process(self.raw)
        measured = measure_region(spectrum, "analyte", self.signal_bounds, self.noise_bounds)
        result = measure_region(spectrum, "analyte", self.signal_bounds, None)
        self.assertEqual(result, replace(measured, noise_bounds_ppm=None, noise_sd=None,
                                         noise_sd_reason="Missing noise window"))
        self.assertGreater(measured.noise_sd, 0)
        self.assertGreater(result.fwhm_hz, 0)
        signal = SignalEvidence(result.name, measurement=result.to_signal_measurement())
        snr = check_snr(self.raw, signal)
        self.assertEqual((snr.status, snr.value, snr.reason), ("cannot be assessed", None, "noise_sd is missing"))
        measured_signal = SignalEvidence(measured.name, measurement=measured.to_signal_measurement())
        self.assertEqual(check_resolution(self.raw, signal), check_resolution(self.raw, measured_signal))

    def test_zero_width_signal_overlaps_noise(self):
        spectrum = process(self.raw)
        result = measure_region(spectrum, "analyte", (2.0, 2.0), self.noise_bounds)
        self.assertIsNone(result.noise_sd)
        self.assertEqual(result.noise_sd_reason, "Noise window overlaps the signal bounds")

    def test_negative_line_keeps_negative_area(self):
        spectrum = process(replace(self.raw, fid=-self.raw.fid))
        result = measure_region(spectrum, "analyte", self.signal_bounds, self.noise_bounds)
        self.assertLess(result.area, 0)
        self.assertGreater(result.noise_sd, 0)
        self.assertIsNone(result.peak_ppm)
        self.assertIsNone(result.height)
        self.assertIsNone(result.fwhm_hz)
        for reason in (result.peak_ppm_reason, result.height_reason, result.fwhm_reason):
            self.assertIn("positive", reason)

    def test_crossing_must_be_inside_bounds_even_between_points(self):
        spectrum = process(self.raw)
        observe = self.raw.observe_frequency_mhz
        narrow = (self.center - 3.9 / observe, self.center + 3.9 / observe)
        result = measure_region(spectrum, "analyte", narrow, self.noise_bounds)
        width = self.rate / np.pi
        self.assertAlmostEqual(result.fwhm_hz, width, delta=0.02 * width)
        for bounds in ((narrow[0], self.center + 3.7 / observe),
                       (self.center - 3.7 / observe, narrow[1])):
            with self.subTest(bounds=bounds):
                result = measure_region(spectrum, "analyte", bounds, self.noise_bounds)
                self.assertIsNone(result.fwhm_hz)
                self.assertIn("crossing", result.fwhm_reason)
                self.assertGreater(result.height, 0)

    def test_line_broadening_blocks_fwhm_only(self):
        spectrum = process(self.raw, ProcessingSettings(lb_hz=2.0))
        result = measure_region(spectrum, "analyte", self.signal_bounds, self.noise_bounds)
        self.assertIsNone(result.fwhm_hz)
        self.assertEqual(result.fwhm_reason, "FWHM requires lb_hz = 0")
        self.assertEqual(result.to_signal_measurement().fwhm_reason, result.fwhm_reason)
        self.assertGreater(result.area, 0)
        self.assertGreater(result.height, 0)
        self.assertGreater(result.noise_sd, 0)

    def test_between_samples_area_does_not_require_a_sampled_peak(self):
        spectrum = process(self.raw)
        bounds = (self.center + 0.0001, self.center + 0.0002)
        result = measure_region(spectrum, "analyte", bounds, self.noise_bounds)
        self.assertGreater(result.area, 0)
        self.assertIsNone(result.peak_ppm)
        self.assertIsNone(result.height)
        self.assertIn("no spectrum points", result.height_reason)

    def test_signal_window_outside_grid_blocks_peak_and_width(self):
        for lb_hz in (0.0, 2.0):
            with self.subTest(lb_hz=lb_hz):
                spectrum = process(self.raw, ProcessingSettings(lb_hz=lb_hz))
                result = measure_region(spectrum, "analyte", (5.4, 12.0), self.noise_bounds)
                for value, reason in ((result.peak_ppm, result.peak_ppm_reason),
                                      (result.height, result.height_reason),
                                      (result.fwhm_hz, result.fwhm_reason)):
                    self.assertIsNone(value)
                    self.assertEqual(reason, "Signal window reaches outside the computed spectral grid")
                self.assertGreater(result.noise_sd, 0)
