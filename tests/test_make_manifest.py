import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]


class MakeManifestTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.deposits_file = self.root / "deposits.json"
        self.out = self.root / "manifest.csv"
        self.deposits = {}

    def add_archive(self, deposit_id, name, files=()):
        deposit = self.deposits.setdefault(deposit_id, {
            "id": deposit_id, "licence": "CC0-1.0", "doi": f"10.0000/{deposit_id}",
            "archives": [],
        })
        folder = self.root / deposit_id / f"x_{Path(name).stem}"
        folder.mkdir(parents=True)
        archive = folder.parent / name
        archive.write_bytes(b"archive")
        deposit["archives"].append({
            "local_name": name, "size_bytes": 7,
            "md5": hashlib.md5(b"archive", usedforsecurity=False).hexdigest(),
        })
        for relative in files:
            path = folder / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(relative.encode("utf-8"))
        return archive

    def run_manifest(self, out=None):
        self.deposits_file.write_text(json.dumps(list(self.deposits.values())), encoding="utf-8")
        return subprocess.run([
            sys.executable, str(PROJECT / "tools/make_manifest.py"),
            "--root", str(self.root), "--deposits", str(self.deposits_file), "--out", str(out or self.out),
        ], capture_output=True, text=True)

    def test_rows_kinds_sorting_and_repeatability(self):
        kinds = {
            "trial/résumé,試.JDF": "jdf",
            "10/fid": "bruker_fid", "10/ser": "bruker_ser",
            "10/acqus": "bruker_acqus", "10/acqu2s": "bruker_acqu2s",
            "10/pdata/1/procs": "bruker_procs", "10/pdata/1/1r": "bruker_1r",
            "sheets/a.xlsx": "workbook", "sheets/b.xls": "workbook",
            "sheets/c.numbers": "workbook", "sheets/d.csv": "workbook",
            "sheets/.DS_Store": "os_junk", "._spectrum.jdf": "os_junk",
            "nested/__MACOSX/10/fid": "os_junk", "read\rme.txt": "other",
        }
        self.add_archive("z", "z.zip", kinds)
        self.add_archive("z", "a.zip", ["10/fid"])
        self.add_archive("a", "raw.zip", ["last.txt"])
        expected = {("z", "z.zip", path): kind for path, kind in kinds.items()}
        expected["z", "a.zip", "10/fid"] = "bruker_fid"
        expected["a", "raw.zip", "last.txt"] = "other"

        result = self.run_manifest()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.count("OK "), 3)
        self.assertIn("z: 16 files", result.stdout)
        self.assertIn("bruker_fid=2", result.stdout)
        self.assertIn("os_junk=3", result.stdout)
        first = self.out.read_bytes()
        self.assertNotIn(b"\r\n", first)
        with self.out.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            self.assertEqual(reader.fieldnames, [
                "deposit_id", "archive", "path", "size_bytes", "sha256", "kind", "licence", "doi",
            ])
            rows = list(reader)
        keys = [(row["deposit_id"], row["archive"], row["path"]) for row in rows]
        self.assertEqual(keys, sorted(expected, key=lambda key: (key[0], key[2], key[1])))
        for row, key in zip(rows, keys):
            payload = row["path"].encode("utf-8")
            self.assertEqual(row["size_bytes"], str(len(payload)))
            self.assertEqual(row["sha256"], hashlib.sha256(payload).hexdigest())
            self.assertEqual(row["kind"], expected[key])
            self.assertEqual(row["licence"], "CC0-1.0")
            self.assertEqual(row["doi"], f"10.0000/{row['deposit_id']}")

        result = self.run_manifest()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.out.read_bytes(), first)

    def test_archive_failures_do_not_write_a_manifest(self):
        archive = self.add_archive("a", "bad.zip")
        self.add_archive("a", "good.zip")
        cases = [("MISMATCH", b"corrupt", 7), ("MISMATCH", b"archive", 5), ("MISSING", None, 7)]
        for status, content, size in cases:
            with self.subTest(content=content):
                self.deposits["a"]["archives"][0]["size_bytes"] = size
                if content is None:
                    archive.unlink()
                else:
                    archive.write_bytes(content)
                result = self.run_manifest()
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual(result.stdout.splitlines()[0].split()[:2], [status, "a/bad.zip"])
                self.assertIn("OK a/good.zip", result.stdout)
                self.assertFalse(self.out.exists())

    def test_malformed_metadata_stops_before_archive_checks(self):
        self.add_archive("a", "good.zip")
        self.add_archive("b", "raw.zip")
        deposit = self.deposits["b"]
        archive = deposit["archives"][0]
        cases = [
            (deposit, "id", "../b"), (deposit, "doi", None),
            (deposit, "licence", ""), (deposit, "licence", "\ud800"),
            (deposit, "archives", []),
            (archive, "local_name", "/raw.zip"), (archive, "md5", "xyz"),
            (archive, "size_bytes", True),
        ]
        for record, field, value in cases:
            with self.subTest(field=field):
                original = record[field]
                record[field] = value
                result = self.run_manifest()
                record[field] = original
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn(field, result.stderr)
                self.assertEqual(result.stdout, "")
                self.assertFalse(self.out.exists())

    def test_output_cannot_overlap_inputs(self):
        archive = self.add_archive("a", "raw.zip", ["fid"])
        folder = archive.parent / "x_raw"
        inside = folder / "manifest.csv"
        alias = self.root / "alias"
        alias.symlink_to(folder, target_is_directory=True)
        for out in (inside, alias / inside.name, archive, self.deposits_file):
            with self.subTest(out=out):
                result = self.run_manifest(out)
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn("--out", result.stderr)
                self.assertEqual(result.stdout, "")
        self.assertFalse(inside.exists())
        self.assertEqual(archive.read_bytes(), b"archive")
        self.assertEqual(json.loads(self.deposits_file.read_text(encoding="utf-8")),
                         list(self.deposits.values()))
