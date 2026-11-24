import copy
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tools.check_case_manifest import validate

PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "data/case_manifest.json"
INVENTORY = PROJECT / "data/manifest.csv"


class CheckCaseManifestTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.manifest = Path(temp.name) / "cases.json"
        self.inventory = Path(temp.name) / "inventory.csv"
        self.data = dict(runs=[], preparations=[], links=[], cases=[], exposures=[], pilot_pairs=[], cohorts={})
        self.processed = ["deposit", "archive", "processed"]
        with self.inventory.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("deposit_id", "archive", "path", "sha256"))
            writer.writerow([*self.processed, "d" * 64])
            for name, deposit in (("a", "deposit"), ("b", "deposit"), ("c", "other")):
                raw = [deposit, "archive", name]
                digest = name * 64
                writer.writerow([*raw, digest])
                self.data["runs"].append(dict(id=name, raw_key=raw, sha256=digest, acquisition_id=None,
                                              alias_of=None, preparation_id=name, processed_key=None,
                                              fixture_aliases=[]))
                self.data["preparations"].append({"id": name})
                self.data["cases"].append(dict(id=name, split="development", run_ids=[name]))
            writer.writerow(["deposit", "archive", "duplicate", "a" * 64])

    def check(self, data=None, inventory=None):
        self.manifest.write_text(json.dumps(self.data if data is None else data), encoding="utf-8")
        inventory = inventory or self.inventory
        return validate(self.manifest, inventory, reference_inventory=inventory)

    def reject(self, message, data=None):
        with self.assertRaisesRegex(ValueError, message):
            self.check(data)

    def test_committed_manifest_and_deterministic_components(self):
        components = validate(MANIFEST, INVENTORY)
        data = json.loads(MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(sum(map(len, components)), 1555 + 86 + 20)
        for section in ("runs", "cases", "preparations", "links"):
            data[section].reverse()
        self.assertEqual(self.check(data, INVENTORY), components)

    def test_each_edge_prevents_split(self):
        edges = {
            "case runs": lambda a, b, data: data["cases"][0]["run_ids"].append("b"),
            "raw sha256": lambda a, b, data: b.update(raw_key=["deposit", "archive", "duplicate"],
                                                       sha256=a["sha256"]),
            "processed to raw": lambda a, b, data: a.update(processed_key=b["raw_key"]),
            "shared processed": lambda a, b, data: (a.update(processed_key=self.processed),
                                                    b.update(processed_key=self.processed)),
            "acquisition": lambda a, b, data: (a.update(acquisition_id={"start": 1, "end": 2}),
                                               b.update(acquisition_id={"end": 2.0, "start": 1})),
            "alias": lambda a, b, data: b.update(alias_of="a"),
            "preparation": lambda a, b, data: b.update(preparation_id="a"),
            "possible link": lambda a, b, data: data["links"].append(
                dict(preparation_ids=["a", "b"], status="possible")),
            "confirmed link": lambda a, b, data: data["links"].append(
                dict(preparation_ids=["a", "b"], status="confirmed")),
            "fixture alias": lambda a, b, data: (a["fixture_aliases"].append("tests/fixtures/example"),
                                                 b["fixture_aliases"].append("tests/fixtures/example")),
        }
        for name, add_edge in edges.items():
            with self.subTest(edge=name):
                data = copy.deepcopy(self.data)
                add_edge(*data["runs"][:2], data)
                self.check(data)
                data["cases"][1]["split"] = "final_holdout"
                self.reject("Split conflict", data)

    def test_excluded_case_keeps_transitive_edges(self):
        self.data["cases"].append(dict(id="bridge", split="excluded", run_ids=["a", "c"]))
        self.data["runs"][2]["alias_of"] = "b"
        self.data["cases"][1]["split"] = "final_holdout"
        self.reject("Split conflict")

    def test_acquisition_is_scoped_to_deposit(self):
        self.data["runs"][0]["acquisition_id"] = self.data["runs"][2]["acquisition_id"] = {"date": 1}
        self.data["cases"][2]["split"] = "final_holdout"
        self.assertIn((("case", "c"), ("preparation", "c"), ("run", "c")), self.check())

    def test_invalid_file_keys_and_references(self):
        for field, value, message in (
            ("raw_key", ["deposit", "archive", "missing"], "raw_key missing"),
            ("sha256", "e" * 64, "raw_key sha256 mismatch"),
            ("processed_key", ["deposit", "archive", "missing"], "processed_key missing"),
            ("processed_sha256", "e" * 64, "processed_key sha256 mismatch"),
            ("alias_of", "missing", "Dangling run reference"),
            ("preparation_id", "missing", "Dangling preparation reference"),
        ):
            with self.subTest(field=field):
                data = copy.deepcopy(self.data)
                data["runs"][0]["processed_key"] = self.processed
                data["runs"][0][field] = value
                self.reject(message, data)
        self.data["links"] = [{"preparation_ids": ["a", "missing"]}]
        self.reject("Dangling preparation reference")
        self.data["links"].clear()
        self.data["cases"][0]["run_ids"] = ["missing"]
        self.reject("Dangling run reference")

    def test_duplicate_ids(self):
        for section in ("runs", "cases", "preparations"):
            with self.subTest(section=section):
                data = copy.deepcopy(self.data)
                data[section].append(data[section][0])
                self.reject("Duplicate .* id", data)

    def test_changed_processed_inventory_digest(self):
        data = json.loads(MANIFEST.read_text(encoding="utf-8"))
        key = next(run["processed_key"] for run in data["runs"] if run["processed_key"] is not None)
        with INVENTORY.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            fields, rows = reader.fieldnames, list(reader)
        row = next(row for row in rows if [row["deposit_id"], row["archive"], row["path"]] == key)
        row["sha256"] = "0" * 64
        with self.inventory.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fields)
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(ValueError, "processed_key sha256 mismatch"):
            validate(MANIFEST, self.inventory)

    def test_fixture_window_requires_owner(self):
        path = "tests/fixtures/example"
        window = self.data["cases"][0]
        window.update(kind="window", path=path)
        self.reject("Fixture window .* has no owner")
        self.data["runs"][1]["fixture_aliases"].append(path)
        self.reject("Fixture window .* has no owner")
        window["run_ids"].append("b")
        self.check()

    def test_pilot_case_requires_ordered_sample_and_ec(self):
        sample, ec = self.data["runs"][:2]
        pair = dict(id="pilot", sample=dict(run_id="a", raw_key=sample["raw_key"]),
                    ec=dict(run_id="b", raw_key=ec["raw_key"]))
        self.data["pilot_pairs"].append(pair)
        case = dict(id="pilot", split="development", run_ids=["a", "b"])
        self.data["cases"].append(case)
        self.check()
        for run_ids in (["a"], ["b", "a"]):
            with self.subTest(run_ids=run_ids):
                case["run_ids"] = run_ids
                self.reject("must match its case's ordered run_ids")
        case["run_ids"] = ["a", "b"]
        del pair["ec"]
        self.reject("Invalid ec run reference")

    def test_unknown_split(self):
        self.data["cases"][0]["split"] = "holdout"
        self.reject("Unknown split")

    def test_holdout_rejects_null_preparation_through_alias(self):
        self.data["cases"][1]["split"] = "excluded"
        self.data["cases"][0]["split"] = "final_holdout"
        self.data["runs"][1].update(alias_of="a", preparation_id=None)
        self.reject("null preparation_id")

    def test_exposure_requires_resolution_and_stays_out_of_holdout(self):
        self.data["cases"][0]["split"] = "final_holdout"
        exposure = {"source": "test", "scope": "worksheet"}
        self.data["exposures"].append(exposure)
        self.reject("unresolved exposure")
        for references, message in (([], "unresolved exposure"), (["missing"], "Dangling exposure"),
                                    (["a"], "exposed run")):
            with self.subTest(references=references):
                exposure["run_ids"] = references
                self.reject(message)
        exposure["run_ids"] = ["b"]
        self.check()

    def test_cli(self):
        command = [sys.executable, str(PROJECT / "tools/check_case_manifest.py")]
        passed = subprocess.run(command, check=True, capture_output=True, text=True)
        self.assertEqual(len(passed.stdout.splitlines()), 1)
        self.assertIn("case manifest: OK", passed.stdout)
        self.manifest.write_text("{", encoding="utf-8")
        failed = subprocess.run([*command, "--manifest", str(self.manifest)], capture_output=True, text=True)
        self.assertNotEqual(failed.returncode, 0)
        self.assertEqual(len(failed.stderr.splitlines()), 1)
        self.assertIn("case manifest: INVALID", failed.stderr)

    def test_non_object_exposure_cli(self):
        self.data["exposures"] = [None]
        self.manifest.write_text(json.dumps(self.data), encoding="utf-8")
        command = [sys.executable, str(PROJECT / "tools/check_case_manifest.py"),
                   "--manifest", str(self.manifest), "--inventory", str(self.inventory)]
        failed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(failed.returncode, 1)
        self.assertEqual(failed.stderr.splitlines(), ["case manifest: INVALID (Exposure must be an object)"])
