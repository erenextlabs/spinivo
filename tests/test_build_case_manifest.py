import csv
import json
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from nmrglue.fileio import jeol

from tools.build_case_manifest import build

PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "data/manifest.csv"
WINDOWS = PROJECT / "data/baseline_windows.json"
JEOL = "tests/fixtures/jeol/nishizaki_fid02_ma_dmso2_d20s_ns32.jdf"
BRUKER = "tests/fixtures/bruker/zenodo17083860_833-1/10"


class BuildCaseManifestTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.temp = Path(temp.name)
        self.root = self.temp / "veri"

    def test_fixture_headers_and_reproducible_cli(self):
        with MANIFEST.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            fields, rows = reader.fieldnames, list(reader)
        jdf = next(row for row in rows if row["sha256"].startswith("dd9aeb9c"))
        fid = next(row for row in rows if row["sha256"].startswith("52fae41b"))
        for row, source, name in ((jdf, JEOL, Path(jdf["path"]).name),
                                  (fid, BRUKER + "/fid", "fid"), (fid, BRUKER + "/acqus", "acqus")):
            folder = self.root / row["deposit_id"] / ("x_" + Path(row["archive"]).stem)
            target = folder / Path(row["path"]).with_name(name)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(PROJECT / source, target)

        alias = dict(jdf, path="duplicate.jdf")
        ser = dict(fid, path=str(Path(fid["path"]).with_name("ser")), kind="bruker_ser")
        orphan = dict(jdf, path="orphan-1-2.jdf")
        wrong_fid = dict(fid, path="earlier/fid")
        wrong_acqus = dict(fid, path="earlier/acqus", kind="bruker_acqus")
        manifest = self.temp / "manifest.csv"
        with manifest.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fields)
            writer.writeheader()
            writer.writerows([alias, ser, orphan, wrong_fid, wrong_acqus, *reversed(rows)])

        result = build(manifest, WINDOWS, self.root)
        runs = {run["raw_key"][2]: run for run in result["runs"]}
        jeol = runs[jdf["path"]]
        self.assertEqual(jeol["processed_key"][2], jdf["path"].replace("-1-1.jdf", "-1-2.jdf"))
        self.assertNotIn(jeol["processed_key"][2], runs)
        self.assertEqual(jeol["acquisition_id"], {
            "start": 929835739, "end": 929836774,
            "title": "MAandDMSO2_T1A_5mm_4226u1534u6000uL_yuzo190619",
        })
        self.assertEqual(jeol["fixture_aliases"], [JEOL])
        self.assertEqual(runs[alias["path"]]["alias_of"], jeol["id"])
        self.assertIsNone(runs[alias["path"]]["acquisition_id"])
        self.assertIsNone(runs[orphan["path"]]["processed_key"])

        bruker = runs[fid["path"]]
        self.assertEqual(bruker["acquisition_id"], {"date": 1690539009})
        self.assertEqual(bruker["processed_key"][2], str(Path(fid["path"]).parent / "pdata/1/1r"))
        self.assertEqual(bruker["fixture_aliases"], [BRUKER])
        self.assertEqual(runs[wrong_fid["path"]]["fixture_aliases"], [])
        self.assertNotIn(ser["path"], runs)

        out = self.temp / "out.json"
        command = [sys.executable, "-E", "-B", str(PROJECT / "tools/build_case_manifest.py"),
                   "--manifest", str(manifest), "--data-root", str(self.root), "--out", str(out)]
        expected = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode("utf-8")
        for _ in range(2):
            subprocess.run(command, check=True, cwd=self.temp, capture_output=True)
            self.assertEqual(out.read_bytes(), expected)

        inside = self.root / "out.json"
        rejected = subprocess.run([*command[:-1], str(inside)], cwd=self.temp, capture_output=True, text=True)
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("--out must be outside --data-root", rejected.stderr)
        self.assertFalse(inside.exists())

        acqus_path = str(Path(fid["path"]).with_name("acqus"))
        manifest.write_text("".join(line for line in manifest.read_text(encoding="utf-8").splitlines(keepends=True)
                                    if acqus_path not in line), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "No matching inventory run for window .*bruker"):
            build(manifest, WINDOWS, self.root)

    def _build_single_row(self, **overrides):
        row = dict(deposit_id="deposit", archive="data.zip", path="run.jdf", kind="jdf", sha256="0" * 64)
        row.update(overrides)
        manifest = self.temp / "manifest.csv"
        with manifest.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, row)
            writer.writeheader()
            writer.writerow(row)
        return build(manifest, WINDOWS, self.root)

    def test_inventory_key_requires_plain_components(self):
        for field, value in (
            ("path", "folder/../run.jdf"), ("path", "/run.jdf"), ("path", "folder//run.jdf"),
            ("path", "./run.jdf"), ("path", r"folder\run.jdf"), ("path", "run.jdf/"),
            ("deposit_id", ".."), ("deposit_id", "folder/deposit"), ("archive", "../data.zip"),
            ("deposit_id", "D:"), ("path", "C:/run.jdf"),
        ):
            with self.subTest(field=field, value=value), self.assertRaisesRegex(
                ValueError, "Invalid inventory key .*" + re.escape(repr(value))
            ):
                self._build_single_row(**{field: value})

    def test_metadata_symlinks_cannot_leave_archive(self):
        folder = self.root / "deposit/x_data"
        folder.mkdir(parents=True)
        outside = folder.parent / "outside"
        outside.write_bytes(b"")
        for kind, raw, metadata in (("jdf", "run.jdf", "run.jdf"), ("bruker_fid", "fid", "acqus")):
            (folder / metadata).symlink_to(outside)
            key = ("deposit", "data.zip", raw)
            with self.subTest(kind=kind), self.assertRaisesRegex(
                ValueError, "Metadata path leaves archive folder for inventory key " + re.escape(repr(key))
            ):
                self._build_single_row(kind=kind, path=raw)

        outside_archive = self.temp / "outside_archive"
        outside_archive.mkdir()
        shutil.copyfile(PROJECT / JEOL, outside_archive / "run.jdf")
        folder.rename(folder.with_name("x_previous"))
        folder.symlink_to(outside_archive, target_is_directory=True)
        key = ("deposit", "data.zip", "run.jdf")
        with self.subTest(kind="archive"), self.assertRaisesRegex(
            ValueError, "Metadata path leaves archive folder for inventory key " + re.escape(repr(key))
        ):
            self._build_single_row()

    def test_jdf_string_scaler_checked_before_parsing(self):
        path = self.root / "deposit/x_data/run.jdf"
        path.parent.mkdir(parents=True)
        data = bytearray((PROJECT / JEOL).read_bytes())
        _, header = jeol.read_header(jeol.IOBuffer(data, jeol.ConversionTable))
        order = "<" if header["endian"] == "little_endian" else ">"
        start = header["param_start"]
        end = start + header["param_length"]
        record = next(offset for offset in range(start + 16, end, 64)
                      if struct.unpack_from(order + "i", data, offset + 32)[0] == 0)
        struct.pack_into(order + "h", data, record + 4, 20)
        path.write_bytes(data)
        with (patch.object(jeol, "read_parameters") as parse,
              self.assertRaisesRegex(ValueError, re.escape(str(path)) + ": invalid JDF parameter scaler")):
            self._build_single_row()
        parse.assert_not_called()

    def test_snapshot_hashes_cases_and_links(self):
        result = build(MANIFEST, WINDOWS, self.root)
        self.assertEqual(len(result["runs"]), 1555)
        self.assertTrue(all(run["acquisition_id"] is None for run in result["runs"]))
        expected = {
            "01": (56, "9bfca77ca06a16ab62026940962503a03758b0d1f3288b773fd2059fcd0f8574"),
            "09": (58, "19503cff5ab8be126deeb5751dc97e6b3eb22523c7a97fc1d198c3e4675f0510"),
            "July": (32, "11f6ef9c493117697469de9ac2b1c030683de1a9a36f064ee6ffd65bea86e8a3"),
        }
        self.assertEqual({name: (len(row["members"]), row["sha256"])
                          for name, row in result["cohorts"].items()}, expected)
        cases = result["cases"]
        windows = json.loads(WINDOWS.read_text(encoding="utf-8"))
        self.assertEqual([{key: case[key] for key in window} for case, window in zip(cases, windows)], windows)
        self.assertEqual(len(cases), 86)
        self.assertTrue(all(case["split"] == "development" for case in cases))
        self.assertEqual(cases[82]["id"], "fid05_widened_dmso2")
        self.assertEqual(cases[82]["bounds"], [2.93, 3.48])
        snapshot = json.loads((PROJECT / "data/case_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(len(result["pilot_pairs"]), 3)
        self.assertEqual(result["pilot_pairs"], snapshot["pilot_pairs"])
        self.assertEqual(len(result["preparations"]), 20)
        self.assertEqual([link["preparation_ids"] for link in result["links"]], [
            ["01:4226u1534u6000uL", "05:4226u1534u2000uL"],
            ["02:3913u2894u6000uL", "06:3913u2894u2000uL"],
            ["03:4280u5692u6000uL", "07:4280u5692u2000uL"],
            ["04:6435u15130u9000uL", "08:6435u15130u3000uL"],
            ["MAandCAF:6435u15130u9000uL", "04:6435u15130u9000uL"],
            ["12:2346u2000uL", "16:2346u2000uL"],
        ])
        self.assertTrue(all(link["status"] == "possible" and link["source"] for link in result["links"]))
        self.assertEqual([exposure["scope"] for exposure in result["exposures"]], [
            "FID03 Sheet1 column E area cells for solutions 09–11",
            "FID02 5 mm sheet rows 146–149, including L147 and L148",
        ])
        for run in snapshot["runs"]:
            run["acquisition_id"] = None
        self.assertEqual(result, snapshot)

    def test_missing_link_preparation_fails_build(self):
        manifest = self.temp / "missing.csv"
        manifest.write_text("".join(line for line in MANIFEST.read_text(encoding="utf-8").splitlines(keepends=True)
                                    if "_4226u1534u2000uL_" not in line), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Missing preparation 05:4226u1534u2000uL"):
            build(manifest, WINDOWS, self.root)

    def test_changed_bound_breaks_cohort_lock(self):
        windows = json.loads(WINDOWS.read_text(encoding="utf-8"))
        windows[0]["bounds"][0] += 0.01
        changed = self.temp / "changed.json"
        changed.write_text(json.dumps(windows), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Cohort 01 hash mismatch"):
            build(MANIFEST, changed, self.root)
