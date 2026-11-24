import unittest
import warnings
from dataclasses import replace
from math import cos, inf, log, nan, nextafter, radians
from pathlib import Path

import numpy as np

from spinivo import RawFID
from spinivo.acquisition import (
    PILOT_THRESHOLDS, RecoveryEvidence, SignalEvidence, SignalMeasurement, ThresholdProfile, check_acquisition,
    check_acquisition_pair, check_recovery, check_resolution, check_snr,
)


class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.raw = RawFID(
            vendor="bruker", source=Path("sample"), nucleus="1H", observe_frequency_mhz=512.0,
            spectral_width_hz=128.0, carrier_ppm=4.6, complex_points=128, acquisition_time_s=1.0,
            relaxation_delay_s=4.0, pulse_width_us=10.0, flip_angle_deg=90.0, pulse_90_us=10.0,
            scans=32, dummy_scans=4, receiver_gain=1.0, temperature_k=298.15, pulse_program="zg",
            digital_filter=None, intensity_note="", fid=np.zeros(128, dtype=np.complex128),
            source_fields={"scans": "NS", "dummy_scans": "DS"}, missing=(),
        )
        self.angle = RecoveryEvidence(calibrated_flip_angle_deg=90.0, calibration_source="pulse log")
        self.width = RecoveryEvidence(calibrated_90_us=10.0, calibration_source="pulse log")
        self.measurement = SignalMeasurement(height=150.0, noise_sd=1.0, fwhm_hz=5.0,
                                             noise_bounds_ppm=(4.5, 4.5625), transformed_points=256)
        self.signal = SignalEvidence("analyte", 1.0, "T1 log", self.measurement)

    def test_thresholds(self):
        for side, status in ((-inf, "flag"), (None, "pass"), (inf, "pass")):
            with self.subTest(side=side):
                tr, height, fwhm = (v if side is None else nextafter(v, side) for v in (log(100), 150., 5.))
                raw = replace(self.raw, relaxation_delay_s=tr - 1.0)
                signal = replace(self.signal, measurement=replace(self.measurement, height=height, fwhm_hz=fwhm))
                recovery, snr, resolution = check_acquisition(raw, [signal], self.angle)
                self.assertEqual([r.status for r in (recovery, snr, resolution)], [status] * 3)
                self.assertAlmostEqual(recovery.value, 0.01)
                self.assertEqual((snr.value, resolution.value), (height, fwhm))

    def test_profile_requires_positive_integer_noise_points(self):
        for value in (0, -1, nan, 32.0, True, None):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "min_noise_points must be a positive integer"):
                ThresholdProfile(min_noise_points=value)
        self.assertEqual(ThresholdProfile(min_noise_points=1).min_noise_points, 1)

    def test_recovery_timing_and_calibration(self):
        short = replace(self.raw, relaxation_delay_s=3.0)
        full = check_recovery(self.raw, self.signal, self.width)
        reduced = check_recovery(short, self.signal, self.width)
        self.assertEqual((full.status, reduced.status), ("pass", "flag"))
        for result in (full, reduced):
            self.assertIn("steady-state bound", result.reason)
            inputs = {i.name: (i.value, i.source) for i in result.inputs}
            self.assertEqual((inputs["scans"], inputs["dummy_scans"]), ((32, "NS"), (4, "DS")))
        self.assertAlmostEqual(full.value, np.exp(-5))
        extra = replace(self.width, extra_recovery_s=1.0, extra_recovery_source="timing log")
        self.assertEqual(check_recovery(short, self.signal, extra).value, full.value)
        thirty = replace(self.width, calibrated_90_us=30.0)
        self.assertEqual(check_recovery(short, self.signal, thirty).status, "pass")
        for program in ("zg30", "zgig30"):
            raw = replace(short, pulse_program=program)
            result = check_recovery(raw, self.signal, self.width)
            self.assertEqual(result.status, "pass")
            inputs = {i.name: (i.value, i.source) for i in result.inputs}
            self.assertEqual(inputs["effective_pulse_width_us"],
                             (10.0 * 0.33, "Bruker zg30/zgig30 pulse program: p1*0.33"))
            self.assertEqual(check_recovery(raw, self.signal, self.angle).value, reduced.value)
            self.assertEqual(check_recovery(replace(raw, vendor="jeol"), self.signal, self.width).value, reduced.value)

    def test_nonfinite_recovery_ratio_is_not_assessed(self):
        for t1 in (1e-320, np.float64(1e-320)):
            with self.subTest(t1_type=type(t1)), np.errstate(over="raise"):
                result = check_recovery(self.raw, replace(self.signal, t1_s=t1), self.angle)
                self.assertEqual((result.status, result.value), ("cannot be assessed", None))
                self.assertEqual(result.reason, "tr_s / t1_s must be finite to evaluate recovery")
        result = check_recovery(replace(self.raw, relaxation_delay_s=1e-320, acquisition_time_s=1e-320),
                                replace(self.signal, t1_s=1e-320), self.angle)
        self.assertEqual(result.status, "flag")
        self.assertAlmostEqual(result.value, np.exp(-2))

    def test_recovery_overflow_is_reported_with_warnings_as_errors(self):
        signal = replace(self.signal, t1_s=np.float64(1e-320))
        with warnings.catch_warnings(), np.errstate(over="warn"):
            warnings.simplefilter("error")
            result = check_recovery(self.raw, signal, self.angle)
        self.assertEqual((result.status, result.value), ("cannot be assessed", None))
        self.assertEqual(result.reason, "tr_s / t1_s must be finite to evaluate recovery")

    def test_zg30_recovery_boundary(self):
        angle = radians(90 * 10.0 * 0.33 / 10.0)
        boundary = log(100 - 99 * cos(angle))
        raw = replace(self.raw, pulse_program="zg30")
        for scale, status in ((1 - 1e-9, "flag"), (1 + 1e-9, "pass")):
            with self.subTest(status=status):
                result = check_recovery(replace(raw, relaxation_delay_s=boundary * scale - 1.0), self.signal, self.width)
                self.assertEqual(result.status, status)
                self.assertAlmostEqual(result.value, 0.01, places=9)

    def test_pulse_program_and_calibration_abstain(self):
        reasons = set()
        for program, vendor in ((None, "bruker"), ("", "bruker"), ("  ", "jeol"),
                                ("cosygpqf", "bruker"), ("proton.jxp", "jeol")):
            with self.subTest(program=program):
                raw = replace(self.raw, vendor=vendor, pulse_program=program)
                result = check_recovery(raw, self.signal, self.angle)
                self.assertEqual((result.status, result.value), ("cannot be assessed", None))
                reasons.add(result.reason)
        self.assertEqual(reasons, {
            "pulse_program is missing, so the excitation scheme is unknown",
            "pulse_program 'cosygpqf' is not a supported 1D pulse-acquire program",
            "pulse_program 'proton.jxp' is not a supported 1D pulse-acquire program",
        })
        for raw in (self.raw, replace(self.raw, vendor="jeol", pulse_program="proton.jxp")):
            result = check_recovery(raw, self.signal)
            self.assertEqual((raw.flip_angle_deg, raw.pulse_90_us, result.status), (90.0, 10.0, "cannot be assessed"))
        self.assertEqual(result.reason, "pulse_program 'proton.jxp' is not a supported 1D pulse-acquire program")
        self.assertEqual(check_recovery(self.raw, self.signal).reason,
                         "calibrated_flip_angle_deg or calibrated_90_us is required")

    def test_required_inputs_and_invalid_values(self):
        for record, evaluate, fields in (
            (self.raw, lambda r: check_recovery(r, self.signal, self.angle),
             {"pulse_program": "ZG", "relaxation_delay_s": 0., "acquisition_time_s": nan, "pulse_width_us": -1.}),
            (self.signal, lambda s: check_recovery(self.raw, s, self.angle), {"t1_s": inf, "t1_source": " "}),
            (self.angle, lambda e: check_recovery(self.raw, self.signal, e),
             {"calibrated_flip_angle_deg": 181., "calibration_source": ""}),
            (self.width, lambda e: check_recovery(self.raw, self.signal, e),
             {"calibrated_90_us": 0., "calibration_source": " "}),
            (self.raw, lambda r: check_resolution(r, self.signal), {"spectral_width_hz": 0., "complex_points": 0}),
            (self.raw, lambda r: check_snr(r, self.signal),
             {"observe_frequency_mhz": inf, "spectral_width_hz": 0., "complex_points": 0}),
            (self.measurement, lambda m: check_resolution(self.raw, replace(self.signal, measurement=m)),
             {"fwhm_hz": 0.}),
            (self.measurement, lambda m: check_snr(self.raw, replace(self.signal, measurement=m)),
             {"height": 0., "noise_sd": -1.}),
        ):
            original = evaluate(record)
            self.assertEqual(original.status, "pass")
            for field, invalid in fields.items():
                for value in (None, invalid):
                    with self.subTest(check=original.check, field=field, value=value):
                        changed = replace(record, **{field: value})
                        result = evaluate(changed)
                        self.assertEqual(result.status, "cannot be assessed")
                        self.assertIsNone(result.value)
                        self.assertIn(field, result.reason)
                        self.assertEqual(evaluate(replace(changed, **{field: getattr(record, field)})), original)

    def test_resolution_reports_missing_fwhm_reason(self):
        for reason, expected in (
            (None, "fwhm_hz is missing"),
            ("FWHM requires lb_hz = 0", "fwhm_hz is missing: FWHM requires lb_hz = 0"),
        ):
            with self.subTest(reason=reason):
                measurement = replace(self.measurement, fwhm_hz=None, fwhm_reason=reason)
                signal = replace(self.signal, measurement=measurement)
                result = check_resolution(self.raw, signal)
                self.assertEqual(result.status, "cannot be assessed")
                self.assertEqual(result.reason, expected)

    def test_noise_window_and_acquired_grid(self):
        before = check_acquisition(self.raw, [self.signal], self.angle)
        for bounds in (None, (4.5, 4.5), (nan, 4.5), (4.5, 4.8)):
            measurement = replace(self.measurement, noise_bounds_ppm=bounds)
            during = check_acquisition(self.raw, [replace(self.signal, measurement=measurement)], self.angle)
            self.assertEqual((during[0], during[2]), (before[0], before[2]))
            self.assertEqual(during[1].status, "cannot be assessed")
            self.assertIn("noise_bounds_ppm", during[1].reason)
            restored = replace(measurement, noise_bounds_ppm=self.measurement.noise_bounds_ppm)
            self.assertEqual(check_snr(self.raw, replace(self.signal, measurement=restored)), before[1])
        raw = replace(self.raw, complex_points=16, acquisition_time_s=0.125, fid=self.raw.fid[:16])
        snr = check_snr(raw, self.signal)
        self.assertEqual(snr.status, "cannot be assessed")
        inputs = {i.name: i.value for i in snr.inputs}
        self.assertEqual((inputs["window_hz"], inputs["noise_points"]), (32.0, 4.0))

    def test_snr_uses_profile_noise_points(self):
        raw = replace(self.raw, complex_points=16, acquisition_time_s=0.125, fid=self.raw.fid[:16])
        profile = replace(PILOT_THRESHOLDS, min_noise_points=4)
        self.assertEqual(check_snr(raw, self.signal, profile).status, "pass")
        result = check_snr(raw, self.signal, replace(profile, min_noise_points=5))
        self.assertEqual(result.status, "cannot be assessed")
        self.assertEqual(result.reason, "noise_points must be at least 5 acquired points")

    def test_pair_preserves_independent_checks(self):
        weak = replace(self.measurement, fwhm_hz=4.0)
        pair = check_acquisition_pair(
            self.raw, [replace(self.signal, t1_s=None, measurement=weak)],
            replace(self.raw, source=Path("calibrant")), [replace(self.signal, name="standard")],
            sample_evidence=self.angle, calibrant_evidence=self.angle,
        )
        self.assertEqual([r.status for r in pair.sample], ["cannot be assessed", "pass", "flag"])
        self.assertEqual([r.status for r in pair.calibrant], ["pass", "pass", "pass"])
        self.assertEqual({(r.run, r.signal) for r in pair.sample}, {("sample", "analyte")})
        self.assertEqual({(r.run, r.signal) for r in pair.calibrant}, {("calibrant", "standard")})
        inputs = {i.name: i.value for i in pair.sample[2].inputs}
        self.assertEqual((inputs["spacing_hz"], inputs["transformed_spacing_hz"]), (1.0, 0.5))
