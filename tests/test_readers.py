import hashlib
import struct
import tempfile
import unittest
from dataclasses import FrozenInstanceError, fields
from pathlib import Path
from unittest.mock import patch

import numpy as np
from nmrglue.fileio import bruker, jeol
from nmrglue.fileio.fileiobase import unit_conversion
from nmrglue.process.proc_base import fft_positive

from spinivo import RawFID, UnsupportedInput, read, read_bruker, read_jeol
from spinivo.acquisition import RecoveryEvidence, SignalEvidence, check_recovery

PROJECT = Path(__file__).resolve().parents[1]
FIXTURES = PROJECT / "tests" / "fixtures"
JEOL = FIXTURES / "jeol" / "nishizaki_fid02_ma_dmso2_d20s_ns32.jdf"
JEOL_ARRAY = FIXTURES / "jeol" / "nishizaki_fid01_pulse_width_array.jdf"
BRUKER = FIXTURES / "bruker" / "zenodo17083860_833-1" / "10"


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))

    def assert_reason(self, reason, path, reader=read):
        with self.assertRaises(UnsupportedInput) as caught:
            reader(path)
        self.assertEqual(caught.exception.reason, reason)
        self.assertIn(str(path), str(caught.exception))
        return str(caught.exception)

    def assert_fid(self, record, expected):
        self.assertIsInstance(record, RawFID)
        self.assertEqual(record.fid.shape, (record.complex_points,))
        self.assertEqual(record.fid.dtype, np.dtype("complex128"))
        self.assertEqual(record.acquisition_time_s, record.complex_points / record.spectral_width_hz)
        np.testing.assert_array_equal(record.fid, expected)

    def assert_fields(self, record, expected):
        for field, value in expected.items():
            with self.subTest(field=field):
                if isinstance(value, float):
                    self.assertAlmostEqual(getattr(record, field), value, places=4)
                else:
                    self.assertEqual(getattr(record, field), value)

    def copy_bruker(self, **files):
        for name, contents in files.items():
            (self.root / name).write_bytes(contents)

    def peak_ppm(self, record):
        obs = record.observe_frequency_mhz
        converter = unit_conversion(
            record.complex_points, True, record.spectral_width_hz, obs, record.carrier_ppm * obs,
        )
        ppm = converter.ppm_scale()
        self.assertTrue(np.all(np.diff(ppm) < 0))
        return ppm[np.argmax(np.abs(fft_positive(record.fid)))]

    def test_fixture_checksums(self):
        # Only the files the reader uses; data/redistribution.json pins acqu and pulseprogram.
        expected = {
            JEOL: "dd9aeb9cafebaa99755ab0576bd8e1c3ac026970766db58e7a777732cd0bcaf3",
            JEOL_ARRAY: "fc9b1dd417bcc15912cacd655cb3891116e5baf26489c05a0e8dac94e9f483f5",
            BRUKER / "acqus": "3aed337ada7a9682461dadf1baab4ee5a3c4dbd55bf221062d285c600c1b11f2",
            BRUKER / "fid": "52fae41b0e075ed323fab8088b87d6a78081f6cb4ed3bdb89cc0171e772b9251",
        }
        for path, digest in expected.items():
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest, str(path))

    def test_jeol_fixture(self):
        record = read(JEOL)
        self.assert_fields(record, {
            "vendor": "jeol", "source": JEOL, "nucleus": "1H",
            "complex_points": 65536, "spectral_width_hz": 7494.0048,
            "observe_frequency_mhz": 399782198.37825 / 1e6,
            "carrier_ppm": 5.0, "relaxation_delay_s": 20.0,
            "pulse_width_us": 7.272, "pulse_90_us": 7.272, "flip_angle_deg": 90.0,
            "scans": 32, "dummy_scans": 4, "receiver_gain": 30,
            "temperature_k": 25.0 + 273.15, "pulse_program": "single_pulse_dec",
            "digital_filter": {"digital_filter": "TRUE", "orders": (2, 54, 73), "factors": (8, 2), "filter_factor": 16},
            "missing": (),
        })
        buffer = JEOL.read_bytes()
        header = jeol.parse_jeol(buffer)["header"]
        self.assertEqual(header["data_format"], "one_d")
        self.assertEqual(header["data_axis_type"][0], "complex")
        self.assertIn("MAandDMSO2", header["title"])
        stored = np.frombuffer(buffer, dtype="<f8", offset=16384, count=2 * 65536)
        self.assert_fid(record, stored[:65536] - 1j * stored[65536:])
        self.assertEqual(record.source_fields["acquisition_time_s"], "x_points, x_sweep")
        self.assertEqual(record.source_fields["pulse_90_us"], "x_90_width")
        self.assertIn("No additional intensity scaling", record.intensity_note)
        # Peak in the matching processed JDF in the deposit.
        self.assertAlmostEqual(self.peak_ppm(record), 2.46256222, delta=0.001)
        with self.assertRaises(FrozenInstanceError):
            record.vendor = "changed"

    def test_bruker_fixture(self):
        record = read_bruker(str(BRUKER))
        self.assert_fields(record, {
            "vendor": "bruker", "source": BRUKER, "nucleus": "1H",
            "complex_points": 32768, "spectral_width_hz": 8417.50841750842,
            "observe_frequency_mhz": 400.33247204,
            "carrier_ppm": 2472.04 / 400.33, "relaxation_delay_s": 30.0,
            # Values from the fixture's acquisition status file.
            "pulse_width_us": 10, "receiver_gain": 42.68, "temperature_k": 298.5974,
            "pulse_90_us": None, "flip_angle_deg": None,
            "scans": 64, "dummy_scans": 4, "pulse_program": "zgig",
            "digital_filter": {"GRPDLY": 67.9872589111328, "DECIM": 2376, "DSPFVS": 20},
            "missing": ("flip_angle_deg", "pulse_90_us"),
        })
        params = bruker.read_jcamp(BRUKER / "acqus")
        self.assertEqual(params["SW"], 21.0262944062838)
        self.assertEqual(params["TD"], 65536)
        self.assertEqual(record.source_fields["acquisition_time_s"], "TD, SW_h")
        self.assertEqual(record.source_fields["carrier_ppm"], "O1 / BF1")
        self.assertEqual(record.source_fields["relaxation_delay_s"], "D[1]")
        self.assertEqual(record.source_fields["pulse_width_us"], "P[1]")
        self.assertIn("NC=-6", record.intensity_note)
        self.assertIn("No NC correction", record.intensity_note)
        stored = np.fromfile(BRUKER / "fid", dtype="<i4")
        self.assert_fid(record, stored[::2] + 1j * stored[1::2])
        # The deposit's procs SF and peaklist.xml give the processed reference and peak.
        reference_shift = (400.330037627532 - 400.33) * 1e6 / record.observe_frequency_mhz
        self.assertAlmostEqual(self.peak_ppm(record) - reference_shift, 1.185527, delta=0.001)

    def test_bruker_distributable_files_only(self):
        # The export excludes acqu and pulseprogram, so the reader must not need them.
        self.copy_bruker(acqus=(BRUKER / "acqus").read_bytes(), fid=(BRUKER / "fid").read_bytes())
        full, minimal = read_bruker(BRUKER), read_bruker(self.root)
        self.assertEqual(minimal.source, self.root)
        np.testing.assert_array_equal(minimal.fid, full.fid)
        for field in fields(RawFID):
            if field.name not in ("source", "fid"):
                with self.subTest(field=field.name):
                    self.assertEqual(getattr(minimal, field.name), getattr(full, field.name))

    def test_jeol_array(self):
        header = jeol.parse_jeol(JEOL_ARRAY.read_bytes())["header"]
        self.assertEqual(header["data_format"], "small_two_d")
        self.assertEqual(header["data_points"][:2], [1024, 12])
        self.assertEqual(header["data_axis_titles"][1], "x_90_width")
        self.assert_reason("arrayed_or_2d", JEOL_ARRAY)

    def test_jeol_processed_spectrum(self):
        processed = self.root / "processed.jdf"
        # First axis type at byte 24, first axis base unit at byte 33.
        for offset, code in ((24, 1), (33, 26), (33, 13)):
            with self.subTest(offset=offset, code=code):
                data = bytearray(JEOL.read_bytes())
                data[offset] = code
                processed.write_bytes(data)
                self.assert_reason("processed_spectrum", processed, read_jeol)

    def test_jeol_numeric_parameters(self):
        path = self.root / "parameter.jdf"
        original = JEOL.read_bytes()
        cases = (
            (b"x_pulse", 0x21, 28, 7.272, "pulse_width_us", 7.272),
            (b"x_pulse", 0x01, 28, 7.272e-6, "pulse_width_us", 7.272),
            (b"X_OFFSET", 0x01, 13, 1998.91099189125, "carrier_ppm", 5.0),
            (b"x_pulse", 0x01, 13, 7.272, None, None),
            (b"x_pulse", 0x21, 28, float("inf"), None, None),
            (b"X_OFFSET", 0x01, 28, 5.0, None, None),
            (b"X_SWEEP", 0x01, 13, 0.0, None, None),
            (b"X_SWEEP", 0x01, 13, -1.0, None, None),
            (b"X_SWEEP", 0x01, 13, float("inf"), None, None),
        )
        for key, prefix_power, unit, value, field, expected in cases:
            with self.subTest(key=key, unit=unit, value=value):
                data = bytearray(original)
                record = original.index(key.ljust(28)) - 36
                data[record + 6:record + 8] = bytes((prefix_power, unit))
                struct.pack_into("<d", data, record + 16, value)
                path.write_bytes(data)
                if field is None:
                    self.assert_reason("truncated", path)
                else:
                    self.assertAlmostEqual(getattr(read(path), field), expected)

    def test_jeol_bad_parameter_records(self):
        path = self.root / "bad-records.jdf"
        original = JEOL.read_bytes()
        for order, offset, value in ((">", 1212, 0), (">", 1216, len(original)),
                                     ("<", 1360, 32), ("<", 1368, 1000), ("<", 1372, 1)):
            with self.subTest(offset=offset):
                data = bytearray(original)
                struct.pack_into(order + "I", data, offset, value)
                path.write_bytes(data)
                with patch.object(jeol, "read_parameters") as parse:
                    self.assert_reason("truncated", path)
                    parse.assert_not_called()

    def test_truncated_jeol(self):
        buffer = JEOL.read_bytes()
        path = self.root / "cut.jdf"
        # Set the float32 type bit while keeping the float64 payload.
        mislabeled = buffer[:14] + bytes([buffer[14] | 0x40]) + buffer[15:]
        for data in (buffer[:8], buffer[:len(buffer) // 2], buffer[:-1], mislabeled):
            path.write_bytes(data)
            self.assert_reason("truncated", path)

    def test_bruker_missing_files(self):
        acqus = (BRUKER / "acqus").read_bytes()
        self.copy_bruker(fid=(BRUKER / "fid").read_bytes())
        self.assert_reason("missing_acqus", self.root)
        self.copy_bruker(acqus=acqus)
        (self.root / "fid").rename(self.root / "ser")
        self.assertIn("without fid", self.assert_reason("ser_not_fid", self.root))
        for evidence, files in (("PARMODE=1", {"acqus": acqus.replace(b"PARMODE= 0", b"PARMODE= 1")}),
                                ("NBL=4", {"acqus": acqus.replace(b"NBL= 1", b"NBL= 4")}),
                                ("acqu2s", {"acqus": acqus, "acqu2s": acqus})):
            with self.subTest(evidence=evidence):
                self.copy_bruker(**files)
                message = self.assert_reason("arrayed_or_2d", self.root)
                self.assertIn(f"found ser; acqus indicates multiple traces or dimensions ({evidence})", message)
        (self.root / "ser").unlink()
        self.assert_reason("missing_fid", self.root)
        (self.root / "acqus").write_bytes(acqus.split(b"##$TE=", 1)[0])
        self.assert_reason("missing_fid", self.root)

    def test_bruker_other_nucleus(self):
        acqus = (BRUKER / "acqus").read_bytes().replace(b"NUC1= <1H>", b"NUC1= <13C>")
        self.copy_bruker(acqus=acqus, fid=(BRUKER / "fid").read_bytes())
        self.assert_reason("not_proton", self.root)

    def test_bruker_nucleus_trailing_comment(self):
        original = (BRUKER / "acqus").read_bytes()
        self.copy_bruker(acqus=original.replace(b"NUC1= <1H>", b"NUC1= <> $$ unassigned"),
                         fid=(BRUKER / "fid").read_bytes())
        self.assertIn("NUC1", self.assert_reason("missing_nucleus", self.root))
        self.copy_bruker(acqus=original.replace(b"NUC1= <1H>", b"NUC1= <1H> $$ observed"))
        self.assertEqual(read(self.root).nucleus, "1H")

    def test_bruker_nucleus_off(self):
        acqus = (BRUKER / "acqus").read_bytes().replace(b"NUC1= <1H>", b"NUC1= <off>")
        self.copy_bruker(acqus=acqus, fid=(BRUKER / "fid").read_bytes())
        self.assertIn("NUC1", self.assert_reason("missing_nucleus", self.root))

    def test_bruker_additional_dimension(self):
        acqus = (BRUKER / "acqus").read_bytes().replace(b"PARMODE= 0", b"PARMODE= 1")
        self.copy_bruker(acqus=acqus, fid=(BRUKER / "fid").read_bytes())
        message = self.assert_reason("arrayed_or_2d", self.root)
        self.assertIn("acqus indicates multiple traces or dimensions (PARMODE=1)", message)
        self.assertNotIn("found ser", message)

    def test_bruker_conflicting_dimension_evidence(self):
        acqus = (BRUKER / "acqus").read_bytes()
        fid = (BRUKER / "fid").read_bytes()
        self.copy_bruker(acqus=acqus, fid=fid, ser=fid, acqu2s=acqus)
        message = self.assert_reason("arrayed_or_2d", self.root)
        self.assertIn("fid and ser", message)
        self.assertIn("acqu2s", message)
        (self.root / "acqu2s").unlink()
        self.assertIn("fid and ser", self.assert_reason("arrayed_or_2d", self.root))
        (self.root / "ser").unlink()
        self.copy_bruker(acqu2s=acqus)
        self.assertIn("acqu2s", self.assert_reason("arrayed_or_2d", self.root))

    def test_bruker_array_value_count(self):
        original = (BRUKER / "acqus").read_bytes()
        fid = (BRUKER / "fid").read_bytes()
        for key, field in (("D", "relaxation_delay_s"), ("P", "pulse_width_us")):
            prefix = f"##${key}= ".encode()
            for bounds in (b"(0..0)", b"(0..40)", b"(0..64)"):
                with self.subTest(key=key, bounds=bounds):
                    self.copy_bruker(acqus=original.replace(prefix + b"(0..63)", prefix + bounds), fid=fid)
                    message = self.assert_reason("truncated", self.root)
                    self.assertTrue(message.endswith(f"parameter {key}"), message)
            start, remainder = original.split(prefix, 1)
            _, tail = remainder.split(b"\n##", 1)
            for array, expected in ((b"(0..0) 7", None), (b"(1..2) 7 8E+0", 7)):
                with self.subTest(key=key, array=array):
                    self.copy_bruker(acqus=start + prefix + array + b"\n##" + tail, fid=fid)
                    self.assertEqual(getattr(read(self.root), field), expected)
            for invalid in (b"bad", b"inf", b"1e309", b"9" * 310):
                with self.subTest(key=key, invalid=invalid):
                    array = b"(0..2) 0 30 " + invalid
                    self.copy_bruker(acqus=start + prefix + array + b"\n##" + tail, fid=fid)
                    message = self.assert_reason("truncated", self.root)
                    self.assertTrue(message.endswith(f"parameter {key}"), message)

    def test_bruker_array_jcamp_number_syntax(self):
        original = (BRUKER / "acqus").read_bytes()
        self.copy_bruker(fid=(BRUKER / "fid").read_bytes())
        for key in ("D", "P"):
            with self.subTest(key=key):
                prefix = f"##${key}= ".encode()
                start, remainder = original.split(prefix, 1)
                _, tail = remainder.split(b"\n##", 1)
                self.copy_bruker(acqus=start + prefix + b"(0..2) 0 30 3_0\n##" + tail)
                self.assertTrue(self.assert_reason("truncated", self.root).endswith(f"parameter {key}"))

    def test_nucleus_diagnoses(self):
        acqus = (BRUKER / "acqus").read_bytes()
        fid = (BRUKER / "fid").read_bytes()
        for value in (b"NUC1= <>", b"NUC1= < >", b"UNUSED= <1H>"):
            with self.subTest(vendor="bruker", value=value):
                self.copy_bruker(acqus=acqus.replace(b"NUC1= <1H>", value), fid=fid)
                self.assertIn("NUC1", self.assert_reason("missing_nucleus", self.root))

        path = self.root / "nucleus.jdf"
        original = JEOL.read_bytes()
        record = original.index(b"X_DOMAIN".ljust(28)) - 36
        for value, key, reason in ((b"", b"X_DOMAIN", "missing_nucleus"),
                                   (b" ", b"X_DOMAIN", "missing_nucleus"),
                                   (b"Proton", b"UNUSED", "missing_nucleus"),
                                   (b"13C", b"X_DOMAIN", "not_proton")):
            with self.subTest(vendor="jeol", value=value, key=key):
                data = bytearray(original)
                data[record + 16:record + 32] = value.ljust(16, b"\x00")
                data[record + 36:record + 64] = key.ljust(28)
                path.write_bytes(data)
                self.assertIn("x_domain", self.assert_reason(reason, path))

    def test_bruker_undecodable_parameter(self):
        acqus = (BRUKER / "acqus").read_bytes()
        fid = (BRUKER / "fid").read_bytes()
        cases = ((b"", b"none", b"<z\x81ig>"), (b"", b"none", b"<z\n\x81ig>"),
                 (b"", "ρ".encode(), b"<z\x81ig>"), (b"\xef\xbb\xbf", "ρ".encode(), b"<z\x81ig>"),
                 (b"", b"\xe9", b"<z\x81ig>"))
        for bom, chemical, pulse_program in cases:
            with self.subTest(bom=bom, chemical=chemical, pulse_program=pulse_program):
                data = acqus.replace(b"CHEMSTR= <none>", b"CHEMSTR= <" + chemical + b">")
                self.copy_bruker(acqus=bom + data.replace(b"PULPROG= <zgig>", b"PULPROG= " + pulse_program),
                                 fid=fid)
                self.assertIn("PULPROG", self.assert_reason("truncated", self.root))

    def test_bruker_undecodable_parameter_after_comment(self):
        acqus = (BRUKER / "acqus").read_bytes()
        acqus = acqus.replace(b"PULPROG= <zgig>", b"PULPROG= <z\n$$ comment\n\x81ig>")
        self.copy_bruker(acqus=acqus, fid=(BRUKER / "fid").read_bytes())
        self.assertIn("PULPROG", self.assert_reason("truncated", self.root))

    def test_bruker_truncated_fid(self):
        self.copy_bruker(acqus=(BRUKER / "acqus").read_bytes(), fid=(BRUKER / "fid").read_bytes()[:-8])
        with patch.object(bruker, "get_data") as decode:
            self.assert_reason("truncated", self.root)
            decode.assert_not_called()

    def test_bruker_truncated_acqus(self):
        original = (BRUKER / "acqus").read_bytes()
        for cut, detail in ((original.replace(b"##$DTYPA= 0", b"##$DTYPA= 1"), "missing or invalid DTYPA"),
                            (original.split(b"##$TE=", 1)[0], "missing JCAMP ##END= line"),
                            (original.split(b"##$D=", 1)[0] + b"##$D= (0..63)\n0",
                             "invalid or incomplete JCAMP parameter D")):
            with self.subTest(detail=detail):
                self.copy_bruker(acqus=cut, fid=(BRUKER / "fid").read_bytes())
                self.assertTrue(self.assert_reason("truncated", self.root).endswith(detail))
        (self.root / "fid").rename(self.root / "ser")
        self.assertTrue(self.assert_reason("truncated", self.root).endswith(detail))

    def test_bruker_malformed_metadata(self):
        original = (BRUKER / "acqus").read_bytes()
        fid = (BRUKER / "fid").read_bytes()
        for key, value in ((b"SW_h", b"bad"), (b"SW_h", b""), (b"SW_h", b"0"),
                           (b"SW_h", b"-1"), (b"SW_h", b"inf"), (b"RG", b"bad"), (b"TE", b"inf"),
                           (b"NS", b"64.0"), (b"PULPROG", b"14")):
            with self.subTest(key=key):
                prefix = b"##$" + key + b"="
                acqus = b"\n".join(prefix + b" " + value if line.startswith(prefix) else line
                                   for line in original.splitlines())
                self.copy_bruker(acqus=acqus, fid=fid)
                self.assert_reason("truncated", self.root)

    def assert_recovery_reason(self, record, expected):
        signal = SignalEvidence("analyte", 1.0, "T1 log")
        result = check_recovery(record, signal, RecoveryEvidence(calibrated_flip_angle_deg=90.0))
        self.assertEqual((result.status, result.value), ("cannot be assessed", None))
        self.assertEqual(result.reason, expected)

    def test_bruker_blank_pulse_program_is_missing(self):
        original = (BRUKER / "acqus").read_bytes()
        for blank in (b"<>", b"< >"):
            with self.subTest(blank=blank):
                self.copy_bruker(acqus=original.replace(b"PULPROG= <zgig>", b"PULPROG= " + blank),
                                 fid=(BRUKER / "fid").read_bytes())
                record = read(self.root)
                self.assertIsNone(record.pulse_program)
                self.assertEqual(record.missing, ("flip_angle_deg", "pulse_90_us", "pulse_program"))
                self.assertNotIn("pulse_program", record.source_fields)
                self.assert_recovery_reason(record, "pulse_program is missing, so the excitation scheme is unknown")

    def test_jeol_blank_and_unenabled_programs(self):
        path = self.root / "program.jdf"
        original = JEOL.read_bytes()
        record = original.index(b"experiment".ljust(28)) - 36
        for value, expected in ((b" " * 16, None), (b"proton.jxp", "proton.jxp")):
            with self.subTest(value=value):
                data = bytearray(original)
                data[record + 16:record + 32] = value.ljust(16, b"\x00")
                path.write_bytes(data)
                raw = read(path)
                self.assertEqual(raw.pulse_program, expected)
                self.assertEqual("pulse_program" in raw.missing, expected is None)
                self.assert_fields(raw, {"flip_angle_deg": 90.0, "pulse_90_us": 7.272, "pulse_width_us": 7.272})
                self.assert_recovery_reason(raw, "pulse_program is missing, so the excitation scheme is unknown"
                                            if expected is None else
                                            "pulse_program 'proton.jxp' is not a supported 1D pulse-acquire program")

    def test_bruker_missing_metadata_and_padded_float_fid(self):
        # Four stored scalars followed by padding to one 1024-byte block.
        scalars = np.zeros(128, dtype=">f8")
        scalars[:4] = [1.5, 2.5, -3.5, 4.5]
        self.copy_bruker(
            acqus=b"##$NUC1= <1H> \n##$AQ_mod= 3\n##$TD= 4\n##$BYTORDA= 1\n##$DTYPA= 2\n##END=\n",
            fid=scalars.tobytes(),
        )
        record = read(self.root)
        self.assertEqual(record.complex_points, 2)
        np.testing.assert_array_equal(record.fid, [1.5 + 2.5j, -3.5 + 4.5j])
        self.assertEqual(set(record.missing), {
            "observe_frequency_mhz", "spectral_width_hz", "carrier_ppm", "acquisition_time_s",
            "relaxation_delay_s", "pulse_width_us", "flip_angle_deg", "pulse_90_us", "scans",
            "dummy_scans", "receiver_gain", "temperature_k", "pulse_program", "digital_filter",
        })
        for field in record.missing:
            self.assertIsNone(getattr(record, field), field)
            self.assertNotIn(field, record.source_fields)
        self.assertIn("NC is absent", record.intensity_note)

    def test_missing_path(self):
        path = self.root / "missing.jdf"
        for reader in (read, read_jeol, read_bruker):
            with self.subTest(reader=reader):
                self.assert_reason("not_found", path, reader)
