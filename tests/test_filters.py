import unittest
from dataclasses import FrozenInstanceError, replace

import numpy as np
from nmrglue.fileio.bruker import bruker_dsp_table

from spinivo import read
from spinivo.filters import remove_delay, resolve_delay
from tests.test_readers import BRUKER, JEOL


class FilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bruker = read(BRUKER)
        cls.jeol = read(JEOL)

    def assert_unresolved(self, raw, field):
        decision = resolve_delay(raw)
        self.assertIsNone(decision.delay_points)
        self.assertEqual(decision.method, "unresolved")
        self.assertIn(field, decision.reason)
        with self.assertRaisesRegex(ValueError, "unresolved digital filter delay"):
            remove_delay(raw, decision)

    def test_bruker_grpdly(self):
        decision = resolve_delay(self.bruker)
        self.assertEqual(decision.delay_points, 67.9872589111328)
        self.assertEqual(decision.method, "GRPDLY")
        self.assertEqual(decision.source_fields, (("GRPDLY", 67.9872589111328),))
        self.assertIsNone(decision.reason)
        with self.assertRaises(FrozenInstanceError):
            decision.delay_points = 0
        raw = replace(self.bruker, digital_filter={"GRPDLY": 3.5, "DSPFVS": 10, "DECIM": 2})
        self.assertEqual(resolve_delay(raw).delay_points, 3.5)

    def test_bruker_legacy_table(self):
        for grpdly in (None, -1, 0, float("inf"), float("nan")):
            with self.subTest(grpdly=grpdly):
                metadata = {"DSPFVS": 10, "DECIM": 2}
                if grpdly is not None:
                    metadata["GRPDLY"] = grpdly
                decision = resolve_delay(replace(self.bruker, digital_filter=metadata))
                self.assertEqual(decision.delay_points, bruker_dsp_table[10][2])
                self.assertEqual(decision.method, "bruker_dsp_table")
                self.assertEqual(dict(decision.source_fields)["DSPFVS"], 10)
                self.assertEqual(dict(decision.source_fields)["DECIM"], 2)

    def test_bruker_unresolved(self):
        for metadata in ({"DSPFVS": 20, "DECIM": 2}, {"GRPDLY": -1, "DSPFVS": 20, "DECIM": 2},
                         {"DSPFVS": 10, "DECIM": 2376}, None):
            with self.subTest(metadata=metadata):
                self.assert_unresolved(replace(self.bruker, digital_filter=metadata), "DSPFVS/DECIM")

    def test_jeol_cascade(self):
        decision = resolve_delay(self.jeol)
        self.assertEqual(decision.delay_points, 19.65625)
        self.assertEqual(decision.method, "jeol_fir_cascade")
        self.assertEqual(dict(decision.source_fields), self.jeol.digital_filter)
        metadata = dict(self.jeol.digital_filter)
        del metadata["filter_factor"]
        raw = replace(self.jeol, digital_filter=metadata)
        self.assertEqual(resolve_delay(raw).delay_points, decision.delay_points)

    def test_jeol_invalid_metadata(self):
        for field, value in (
            ("digital_filter", None), ("digital_filter", "unknown"),
            ("orders", None), ("orders", (3, 54, 73)), ("orders", (2, 54, 0)),
            ("factors", None), ("factors", (8,)), ("factors", (8, 1.5)), ("factors", (8, True)),
            ("filter_factor", 15), ("filter_factor", 16.5), ("filter_factor", None),
        ):
            with self.subTest(field=field, value=value):
                metadata = dict(self.jeol.digital_filter)
                if value is None and field == "digital_filter":
                    del metadata[field]
                else:
                    metadata[field] = value
                self.assert_unresolved(replace(self.jeol, digital_filter=metadata), field)

    def test_jeol_filter_disabled(self):
        raw = replace(self.jeol, digital_filter={"digital_filter": "FALSE"})
        decision = resolve_delay(raw)
        self.assertEqual(decision.delay_points, 0)
        self.assertEqual(decision.method, "filter_disabled")
        correction = remove_delay(raw, decision)
        np.testing.assert_array_equal(correction.fid, raw.fid)
        self.assertEqual(correction.method, "bypass")
        self.assertEqual(correction.tail_trim_points, 0)

    def test_synthetic_without_filter(self):
        raw = replace(self.jeol, vendor="synthetic", digital_filter=None)
        decision = resolve_delay(raw)
        self.assertEqual(decision.delay_points, 0.0)
        self.assertEqual(decision.method, "synthetic_no_filter")
        self.assertIsNone(decision.reason)

    def test_synthetic_with_filter_is_unresolved(self):
        for metadata in ({}, {"digital_filter": False}, {"GRPDLY": 1.0}):
            with self.subTest(metadata=metadata):
                self.assert_unresolved(replace(self.jeol, vendor="synthetic", digital_filter=metadata), "synthetic")

    def test_real_fid_onset(self):
        for raw, tail_trim in ((self.bruker, 68), (self.jeol, 20)):
            with self.subTest(vendor=raw.vendor):
                original = raw.fid.copy()
                before = np.abs(raw.fid[:64])
                decision = resolve_delay(raw)
                correction = remove_delay(raw, decision)
                after = np.abs(correction.fid[:64])
                self.assertLess(before[0], 0.5 * max(before))
                self.assertGreaterEqual(after[0], 0.5 * max(after))
                self.assertEqual(correction.decision, decision)
                self.assertEqual(correction.method, "numpy.fft.fft / numpy.fft.ifft")
                self.assertEqual(correction.first_point_treatment, "unweighted after shift")
                self.assertEqual(correction.tail_trim_points, tail_trim)
                self.assertEqual(len(correction.fid), raw.complex_points - tail_trim)
                np.testing.assert_array_equal(raw.fid, original)

    def test_delay_leaving_no_samples_is_refused(self):
        raw = replace(self.bruker, fid=self.bruker.fid[:8], complex_points=8)
        with self.assertRaisesRegex(ValueError, "no usable FID points"):
            remove_delay(raw, resolve_delay(raw))

    def test_fractional_delay_restores_phase_across_carrier(self):
        size = 4096
        bins = np.array([-1280, -512, 704, 1280])
        phases = np.array([0.2, -0.5, 0.8, -1.1])
        delay = 19.65625
        time = np.arange(size) - delay
        # Sample four damped lines starting at the known fractional onset.
        rates = -1 / 180 + 2j * np.pi * bins / size
        delayed = np.exp(rates[:, None] * time + 1j * phases[:, None]).sum(axis=0)
        delayed[time < 0] = 0
        raw = replace(self.bruker, fid=delayed, complex_points=size, digital_filter={"GRPDLY": delay})
        decision = resolve_delay(raw)

        def phase_errors(fid):
            spectrum = np.fft.fft(fid, n=size)
            peaks = []
            for index in bins % size:
                window = np.abs(spectrum[index - 4:index + 5])
                peaks.append(index - 4 + np.argmax(window))
            return np.abs(np.angle(spectrum[peaks] * np.exp(-1j * phases), deg=True))

        corrected = remove_delay(raw, decision)
        np.testing.assert_array_less(phase_errors(corrected.fid), 2)
        wrong_sign = np.fft.ifft(np.fft.fft(delayed) * np.exp(-2j * np.pi * np.fft.fftfreq(size) * delay))
        integer = remove_delay(raw, replace(decision, delay_points=round(delay)))
        for label, fid in (("uncorrected", delayed), ("wrong sign", wrong_sign), ("integer", integer.fid)):
            with self.subTest(control=label):
                self.assertGreater(max(phase_errors(fid)), 5)
