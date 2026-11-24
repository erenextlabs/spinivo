import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
FIXTURES = PROJECT / "tests" / "fixtures"
ALLOWLIST = PROJECT / "data" / "redistribution.json"
TOOL = PROJECT / "tools" / "export_check.py"
BRUKER = "tests/fixtures/bruker/zenodo17083860_833-1/10"
CREDITS = ["README.md", "tests/fixtures/README.md"]


def repository_files():
    """Allowlisted root and data paths plus the source folders, without caches or untracked root files."""
    rules = json.loads(ALLOWLIST.read_text(encoding="utf-8"))["rules"]
    names = {rule["path"] for rule in rules if "path" in rule and (PROJECT / rule["path"]).is_file()}
    for folder in ("src", "tests", "tools", "docs"):
        names.update(path.relative_to(PROJECT).as_posix() for path in (PROJECT / folder).rglob("*")
                     if path.is_file() and "__pycache__" not in path.parts)
    return sorted(names)


class ExportCheckTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))

    def run_check(self, names, root=PROJECT, null=True):
        separator = "\0" if null else "\n"
        return subprocess.run(
            [sys.executable, TOOL, "--root", root, "--allowlist", ALLOWLIST, "-z" if null else "--stdin"],
            input=separator.join(names).encode(), capture_output=True, check=False)

    def rows(self, result):
        lines = result.stdout.decode().splitlines()
        return {fields[3]: fields for fields in (line.split("\t") for line in lines) if len(fields) == 5}

    def copy(self, names, root=None):
        root = root or self.root
        for name in names:
            (root / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(PROJECT / name, root / name)
        return list(names)

    def test_fixtures_listed_with_hashes(self):
        rules = {rule["path"]: rule for rule in json.loads(ALLOWLIST.read_text(encoding="utf-8"))["rules"]
                 if "path" in rule}
        fixtures = [path for path in FIXTURES.rglob("*") if path.is_file()]
        self.assertGreaterEqual(len(fixtures), 13)
        for path in fixtures:
            name = path.relative_to(PROJECT).as_posix()
            with self.subTest(path=name):
                self.assertIn(name, rules)
                self.assertEqual(rules[name]["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())

    def test_repository_dry_run(self):
        result = self.run_check(repository_files())
        self.assertEqual(result.returncode, 0, result.stdout.decode() + result.stderr.decode())
        rows = self.rows(result)
        self.assertNotIn("unknown", {row[0] for row in rows.values()})
        for name in ("acqu", "pulseprogram"):
            if (PROJECT / BRUKER / name).exists():
                self.assertEqual(rows[f"{BRUKER}/{name}"][:2], ["excluded", "CC-BY-4.0"])
        self.assertEqual(rows[f"{BRUKER}/fid"][:2], ["cleared", "CC-BY-4.0"])
        self.assertEqual(rows["src/spinivo/cli.py"][:2], ["cleared", "GPL-3.0-or-later"])
        self.assertEqual(rows["LICENSE"][:2], ["cleared", "LicenseRef-FSF-GPL-3.0-text (verbatim copying only)"])
        self.assertIn("CC-BY-4.0", rows["tests/fixtures/processing_references.json"][1])

    def test_changed_byte_and_unlisted_file(self):
        names = self.copy([path.relative_to(PROJECT).as_posix() for path in FIXTURES.rglob("*")
                           if path.is_file()] + ["README.md"])
        self.assertEqual(self.run_check(names, root=self.root).returncode, 0)

        fid = self.root / BRUKER / "fid"
        data = bytearray(fid.read_bytes())
        data[100] ^= 1
        fid.write_bytes(data)
        result = self.run_check(names, root=self.root)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(self.rows(result)[f"{BRUKER}/fid"][0], "mismatch")
        shutil.copyfile(PROJECT / BRUKER / "fid", fid)

        for extra in ("tests/fixtures/extra.bin", "tests/fixtures/extra.py"):
            with self.subTest(extra=extra):
                (self.root / extra).write_bytes(b"extra")
                result = self.run_check([*names, extra], root=self.root)
                self.assertEqual(result.returncode, 1)
                self.assertEqual(self.rows(result)[extra][:2], ["unknown", "-"])

        result = self.run_check(["tests/fixtures/absent.json", "../outside"], root=self.root)
        self.assertEqual(result.returncode, 1)
        self.assertEqual({row[0] for row in self.rows(result).values()}, {"missing"})

    def test_symlinks_are_unknown(self):
        names = self.copy(["README.md", "src/spinivo/cli.py"])
        os.symlink(self.root / "README.md", self.root / "docs.md")
        (self.root / "docs").symlink_to(self.root / "src", target_is_directory=True)
        result = self.run_check([*names, "docs.md", "docs/spinivo/cli.py"], root=self.root)
        self.assertEqual(result.returncode, 1)
        rows = self.rows(result)
        self.assertEqual(rows["src/spinivo/cli.py"][0], "cleared")
        for name in ("docs.md", "docs/spinivo/cli.py"):
            self.assertEqual(rows[name][0::4], ["unknown", "symlink"])

    def test_path_list_input(self):
        self.copy(["tools/export_check.py"])
        (self.root / "src").mkdir()
        (self.root / "src/cli.py ").write_bytes(b"spaced")
        for null, names in ((False, ["src/cli.py ", "tools/export_check.py\r", ""]),
                            (True, ["src/cli.py ", "tools/export_check.py", " tools/export_check.py"])):
            with self.subTest(null=null):
                rows = self.rows(self.run_check(names, root=self.root, null=null))
                self.assertEqual(rows["src/cli.py "][2], hashlib.sha256(b"spaced").hexdigest())
                self.assertEqual(rows["tools/export_check.py"][0], "cleared")
                self.assertEqual(len(rows), 3 if null else 2)
        self.assertEqual(rows[" tools/export_check.py"][0], "missing")

    def test_attribution_files_required(self):
        names = self.copy([f"{BRUKER}/acqus", f"{BRUKER}/fid", *CREDITS])
        self.assertEqual(self.run_check(names, root=self.root).returncode, 0)
        result = self.run_check(names[:3], root=self.root)
        self.assertEqual(result.returncode, 1)
        row = self.rows(result)["tests/fixtures/README.md"]
        self.assertEqual(row[0], "attribution")
        self.assertIn(f"attribution missing for {BRUKER}/acqus", row[4])
        jeol = self.copy(["tests/fixtures/jeol/nishizaki_fid02_ma_dmso2_d20s_ns32.jdf"])
        self.assertEqual(self.run_check(jeol, root=self.root).returncode, 0)

    def test_clean_export_runs_reader_tests(self):
        rows = self.rows(self.run_check(repository_files()))
        cleared = [name for name, row in rows.items() if row[0] == "cleared"]
        self.copy(cleared)
        self.assertFalse((self.root / BRUKER / "pulseprogram").exists())
        self.assertEqual(self.run_check(cleared, root=self.root).returncode, 0)
        result = subprocess.run([sys.executable, "-m", "unittest", "tests.test_readers"], cwd=self.root,
                                env={**os.environ, "PYTHONPATH": str(self.root / "src")},
                                capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr[-2000:])


if __name__ == "__main__":
    unittest.main()
