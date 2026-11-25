import io
import json
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest.mock import patch

from tools import check_case_policy, compare_references
from tools.check_case_policy import MANIFEST, PROJECT, WINDOWS, check_cohorts, preflight_purity


class CasePolicyTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.manifest = Path(temp.name) / "cases.json"
        self.data = self.load_manifest()
        self.pairs = [(pair["sample"]["raw_key"], pair["ec"]["raw_key"]) for pair in self.data["pilot_pairs"]]

    def load_manifest(self):
        return json.loads(MANIFEST.read_text(encoding="utf-8"))

    def write_manifest(self, data):
        self.manifest.write_text(json.dumps(data), encoding="utf-8")

    def test_pinned_keys_and_pairs(self):
        for sample, ec in self.pairs:
            preflight_purity("internal_standard", sample)
            preflight_purity("internal_standard", tuple(ec))
            preflight_purity("external_calibration", sample, ec)

    def test_preflight_without_site_packages(self):
        code = ("from tools.check_case_policy import preflight_purity\n"
                f"preflight_purity('external_calibration', *{self.pairs[0]!r})")
        subprocess.run([sys.executable, "-S", "-c", code], cwd=PROJECT, check=True)

    def test_processed_sibling_refused(self):
        sample, ec = self.pairs[0]
        run = next(run for run in self.data["runs"] if run["raw_key"] == sample)
        sibling = run["processed_key"]
        self.assertTrue(sibling[2].endswith("-1-2.jdf"))
        with self.assertRaisesRegex(ValueError, "sample key is not a pinned pilot key"):
            preflight_purity("external_calibration", sibling, ec)

    def test_other_file_in_same_folder_refused(self):
        sample, ec = self.pairs[0]
        folder = Path(sample[2]).parent
        other = next(key for key in (run["raw_key"] for run in self.data["runs"])
                     if key[:2] == sample[:2] and Path(key[2]).parent == folder and key not in (sample, ec))
        with self.assertRaisesRegex(ValueError, "EC key is not a pinned pilot key"):
            preflight_purity("external_calibration", sample, other)

    def test_cross_pair_and_reversed_pair_refused(self):
        sample, ec = self.pairs[0]
        for keys in ((sample, self.pairs[1][1]), (ec, sample)):
            with self.subTest(keys=keys), self.assertRaisesRegex(ValueError, "sample, EC order"):
                preflight_purity("external_calibration", *keys)

    def test_final_holdout_refused(self):
        self.data["exposures"] = []
        for case in self.data["cases"]:
            case["split"] = "final_holdout" if case["id"] == "pilot-FID02" else "excluded"
        self.write_manifest(self.data)
        with (patch.object(check_case_policy, "MANIFEST", self.manifest),
              self.assertRaisesRegex(ValueError, "final_holdout component")):
            preflight_purity("internal_standard", self.pairs[0][0])

    def test_exposed_run_refused(self):
        run_id = self.data["pilot_pairs"][0]["ec"]["run_id"]
        self.data["exposures"].append({"scope": "test", "run_ids": [run_id]})
        self.write_manifest(self.data)
        with (patch.object(check_case_policy, "MANIFEST", self.manifest),
              self.assertRaisesRegex(ValueError, "exposed component")):
            preflight_purity("external_calibration", *self.pairs[0])

    def test_preflight_validates_before_checking_keys(self):
        self.data["exposures"].append(None)
        self.write_manifest(self.data)
        with (patch.object(check_case_policy, "MANIFEST", self.manifest),
              self.assertRaisesRegex(ValueError, "Exposure must be an object")):
            preflight_purity("internal_standard", ("invalid", "key", "here"))

    def test_committed_cohorts(self):
        cohorts = check_cohorts()
        self.assertEqual({name: len(cohort["members"]) for name, cohort in cohorts.items()},
                         {"01": 56, "09": 58, "July": 32})
        self.data["cases"].reverse()
        self.write_manifest(self.data)
        self.assertEqual(check_cohorts(self.manifest), cohorts)

    def test_same_count_window_field_changes_refused(self):
        original = next(case for case in self.data["cases"] if case["kind"] == "window")
        for field, value in (
            ("path", "other.jdf"), ("bounds", [6.16, 6.28]),
            ("digests", dict.fromkeys(original["digests"], "0" * 64)),
            ("settings", dict(original["settings"], lb_hz=1)),
        ):
            with self.subTest(field=field):
                data = self.load_manifest()
                window = next(case for case in data["cases"] if case["kind"] == "window")
                window[field] = value
                self.write_manifest(data)
                with self.assertRaisesRegex(ValueError, "Cohort .* hash mismatch"):
                    check_cohorts(self.manifest)

    def test_same_count_cohort_membership_changes_refused(self):
        members = self.data["cohorts"]["01"]["members"]
        members[-1] = next(member for member in self.data["cohorts"]["09"]["members"] if member not in members)
        self.write_manifest(self.data)
        with self.assertRaisesRegex(ValueError, "cohort records differ"):
            check_cohorts(self.manifest)
        data = self.load_manifest()
        next(case for case in data["cases"] if case["id"] == members[0])["cohorts"] = []
        self.write_manifest(data)
        with self.assertRaisesRegex(ValueError, "cohort membership"):
            check_cohorts(self.manifest)

    def test_widened_case_stays_separate(self):
        for field, value in (("kind", "window"), ("cohorts", ["09"])):
            with self.subTest(field=field):
                data = self.load_manifest()
                widened = next(case for case in data["cases"] if case["id"] == "fid05_widened_dmso2")
                widened[field] = value
                self.write_manifest(data)
                with self.assertRaisesRegex(ValueError, "must be an overlap case outside"):
                    check_cohorts(self.manifest)

    def test_policy_cli(self):
        command = [sys.executable, "-S", str(PROJECT / "tools/check_case_policy.py")]
        passed = subprocess.run(command, check=True, capture_output=True, text=True)
        self.assertEqual(passed.stdout.splitlines(), ["case policy: OK (01: 56, 09: 58, July: 32)"])
        snapshot = json.loads(WINDOWS.read_text(encoding="utf-8"))
        snapshot[0]["bounds"][0] += 0.01
        windows = self.manifest.with_name("windows.json")
        windows.write_text(json.dumps(snapshot), encoding="utf-8")
        failed = subprocess.run([*command, "--windows", str(windows)], capture_output=True, text=True)
        self.assertEqual(failed.returncode, 1)
        self.assertEqual(len(failed.stderr.splitlines()), 1)
        self.assertIn("case policy: INVALID (Cohort 01 hash mismatch", failed.stderr)

        snapshot[0]["path"] = None
        windows.write_text(json.dumps(snapshot), encoding="utf-8")
        failed = subprocess.run([*command, "--windows", str(windows)], capture_output=True, text=True)
        self.assertEqual(failed.returncode, 1)
        self.assertEqual(failed.stderr.splitlines(), ["case policy: INVALID (Baseline window path must be a string)"])

    def test_policy_cli_refuses_missing_or_changed_pair(self):
        command = [sys.executable, str(PROJECT / "tools/check_case_policy.py"), "--manifest", str(self.manifest)]
        for missing in (True, False):
            with self.subTest(missing=missing):
                data = self.load_manifest()
                if missing:
                    data["pilot_pairs"].pop()
                else:
                    pair = data["pilot_pairs"][0]
                    pair["sample"], pair["ec"] = pair["ec"], pair["sample"]
                    next(case for case in data["cases"] if case["id"] == pair["id"])["run_ids"].reverse()
                self.write_manifest(data)
                failed = subprocess.run(command, capture_output=True, text=True)
                self.assertEqual(failed.returncode, 1)
                self.assertEqual(failed.stderr.splitlines(), [
                    "case policy: INVALID (Manifest pilot pairs differ from the pinned pilot pairs)",
                ])

    def test_policy_cli_refuses_swapped_pilot_pairs(self):
        pairs = {pair["id"]: pair for pair in self.data["pilot_pairs"]}
        left, right = pairs["pilot-FID02"], pairs["pilot-FID05"]
        for role in ("sample", "ec"):
            left[role], right[role] = right[role], left[role]
        for case in self.data["cases"]:
            if case["id"] in (left["id"], right["id"]):
                pair = pairs[case["id"]]
                case["run_ids"] = [pair[role]["run_id"] for role in ("sample", "ec")]
        self.write_manifest(self.data)
        command = [sys.executable, str(PROJECT / "tools/check_case_policy.py"), "--manifest", str(self.manifest)]
        failed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(failed.returncode, 1)
        self.assertEqual(failed.stderr.splitlines(), [
            "case policy: INVALID (Manifest pilot pairs differ from the pinned pilot pairs)",
        ])

    def test_compare_reports_all_pilots_in_order(self):
        with (patch.object(sys, "argv", ["compare_references.py", str(self.manifest.parent)]),
              patch.object(compare_references, "report_pilot") as report):
            compare_references.main()
        self.assertEqual([call.args[1]["name"] for call in report.call_args_list], ["FID02", "FID05", "FID06"])

    def test_compare_refuses_missing_or_disagreeing_cross_run_ratios(self):
        for ratios, detail in (
            ([["G5", "H5"]], "Pilot requires a cross-run ratio"),
            ([["G5", "F4"], ["F4", "H5"]], "Cross-run ratios disagree on sample and EC roles"),
        ):
            references = compare_references.load_references()
            references["pilots"][-1]["ratios"] = ratios
            err = io.StringIO()
            with (self.subTest(ratios=ratios),
                  patch.object(sys, "argv", ["compare_references.py", str(self.manifest.parent)]),
                  patch.object(compare_references, "load_references", return_value=references),
                  patch.object(compare_references, "report_pilot") as report,
                  redirect_stderr(err), self.assertRaises(SystemExit) as stopped):
                compare_references.main()
            self.assertEqual(stopped.exception.code, 1)
            report.assert_not_called()
            self.assertEqual(err.getvalue().splitlines(), [f"Pilot FID06 refused: {detail}"])

    def test_compare_refuses_before_reading_any_pilot(self):
        references = compare_references.load_references()
        run = references["pilots"][-1]["runs"][0]
        run["path"] = run["path"].removesuffix("-1-1.jdf") + "-1-2.jdf"
        err = io.StringIO()
        with (patch.object(sys, "argv", ["compare_references.py", str(self.manifest.parent)]),
              patch.object(compare_references, "load_references", return_value=references),
              patch.object(compare_references, "report_pilot") as report,
              redirect_stderr(err), self.assertRaises(SystemExit) as stopped):
            compare_references.main()
        self.assertEqual(stopped.exception.code, 1)
        report.assert_not_called()
        self.assertIn("Pilot FID06 refused: EC key is not a pinned pilot key", err.getvalue())
