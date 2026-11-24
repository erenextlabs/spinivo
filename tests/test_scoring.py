import json
import unittest
from collections import defaultdict
from dataclasses import asdict, replace
from hashlib import sha256
from pathlib import Path
from uuid import uuid4

from spinivo import read
from spinivo.integral_components import ComponentEvidence, Line, check_region
from spinivo.processing import process
from spinivo.scoring import Prediction, evaluate_integrals, score_fixture, score_integrals
from tests.integral_fixtures import (
    FIXTURES, SETTINGS, declared_evidence, load_cases, make_raw, ppm_bounds, processing_regions, target_lines,
)
from tools.compare_references import load_references


PROCESSING_VERSION = "processing-5bf8bfc91ab6"


def effective_lines(lines):
    return tuple(Line(ppm_bounds((center,))[0], area, width + SETTINGS.lb_hz) for center, area, width in lines)


def prepare_fixture(rows):
    """Build traces and declared checker inputs before calling the evaluator."""
    prepared, real_runs, settings_used = {}, defaultdict(list), set()
    synthetic_traces = {}
    by_id = {case["id"]: case for case in rows}
    for case in rows:
        if case["kind"] == "real":
            real_runs[case["source"]].append(case)
            continue
        parent_id = case.get("parent", case["id"])
        parent = by_id[parent_id]
        if parent_id not in synthetic_traces:
            raw = replace(make_raw(parent), source=Path(f"synthetic-{uuid4().hex}"))
            synthetic_traces[parent_id] = raw, process(raw, SETTINGS, regions_ppm=processing_regions(parent))
            settings_used.add(SETTINGS)
        raw, spectrum = synthetic_traces[parent_id]
        satellites = parent.get("target") == "satellites"
        target = target_lines(parent)
        declared = declared_evidence(case)
        foreign = declared.get("extra_line_assignment", ())
        noise = declared.get("noise_window_hz")
        components = ComponentEvidence(
            target_lines=effective_lines(target[1:2] if satellites else target), target_source="synthetic target model",
            foreign_lines=effective_lines(foreign), foreign_source="declared extra-line assignment" if foreign else None,
            satellite_convention="include" if satellites else "absent",
            coupling_hz=130 if satellites else None, satellite_source="synthetic target model",
        )
        prepared[case["id"]] = dict(raw=raw, spectrum=spectrum, bounds=ppm_bounds(case["bounds"]),
                                    components=components, noise_bounds_ppm=ppm_bounds(noise) if noise else None)

    noise_windows = load_references()["noise_windows_ppm"]
    for source, cases in real_runs.items():
        path = FIXTURES.parents[1] / source
        if not path.exists():
            continue
        first = cases[0]
        settings = replace(SETTINGS, reference_shift_ppm=first.get("reference_shift_ppm", 0),
                           reference_source=first.get("reference_source"))
        regions = tuple(tuple(case["bounds"]) for case in cases)
        raw = read(path)
        spectrum = process(raw, settings, regions_ppm=regions)
        settings_used.add(settings)
        for case, bounds in zip(cases, regions):
            prepared[case["id"]] = dict(raw=raw, spectrum=spectrum, bounds=bounds, components=ComponentEvidence(),
                                        noise_bounds_ppm=tuple(noise_windows[raw.vendor]))

    settings_json = sorted(json.dumps(asdict(settings), sort_keys=True) for settings in settings_used)
    digest = sha256("\n".join(settings_json).encode()).hexdigest()[:12]
    return prepared, f"processing-{digest}"


def fixture_predictions(inputs):
    return {
        case["key"]: Prediction.from_region(check_region(
            case["spectrum"], case["bounds"], case["components"], raw=case["raw"], noise_bounds_ppm=case["noise_bounds_ppm"]
        ))
        for case in inputs
    }


def row(name, label, *causes, kind="real", source=None, **metadata):
    return {"id": name, "kind": kind, "source": source or name, "bounds": [1, 2],
            "label": label, "causes": list(causes), **metadata}


class ScoringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = score_fixture(FIXTURES / "integral_cases.json", fixture_predictions,
                                    prepare=prepare_fixture)

    def assert_counts(self, cohort, *, labels, statuses, flags, errors, scores, assessed, groups):
        labelled, determinate = sum(labels), sum(labels[:2])
        self.assertEqual(cohort["labelled"], labelled)
        self.assertEqual(cohort["determinate"], determinate)
        self.assertEqual(cohort["groups"], groups)
        for name, values, keys in (
            ("labels", labels, ("reliable", "unreliable", "indeterminate")),
            ("statuses", statuses, ("flag", "pass", "indeterminate", "not assessed")),
            ("region_flags", flags, ("reliable", "unreliable", "indeterminate")),
        ):
            self.assertEqual(cohort[name], dict(zip(keys, values)), name)
        for keys, values in (
            (("true_positive", "false_positive", "false_negative"), errors),
            (("precision", "recall", "coverage", "determinate_coverage"), scores),
            (("assessed", "determinate_assessed"), assessed),
        ):
            self.assertEqual(tuple(cohort[key] for key in keys), values, keys)

    def assert_causes(self, cohort, expected):
        self.assertEqual(cohort["causes"], [
            {"cause": cause, "expected": known, "flagged": flagged, "matched": matches,
             "flagged_indeterminate": uncertain, "processing_version": version}
            for cause, known, flagged, matches, uncertain, version in expected
        ])

    def test_hand_counted_cohorts_and_source_groups(self):
        cases = [
            row("hit", "unreliable", "bounds", "phase", source="shared", run="first"),
            row("wrong", "unreliable", "baseline", source="shared", run="second"),
            row("unclear", "unreliable", "satellites"),
            row("missing", "unreliable", "baseline"),
            row("false_flag", "reliable"),
            row("pass", "reliable"),
            row("unknown", "indeterminate", "neighbour"),
            row("control", "unreliable", "bounds", source="shared", processing_only=True),
            row("clean", "reliable", kind="synthetic", source="generator"),
            row("weak", "unreliable", "low_snr", kind="synthetic", source="generator", snr=30),
        ]
        variant = row("variant", "unreliable", "satellites", kind="synthetic", source="generator", parent="clean")
        predictions = {
            "hit": Prediction("flag", ("bounds", "phase", "bounds")),
            "wrong": Prediction("flag", ("snr",)),
            "unclear": Prediction("cannot be assessed"),
            "false_flag": Prediction("flag", ("baseline",)),
            "pass": Prediction("pass"),
            "unknown": Prediction("flag", ("neighbour",)),
            "control": Prediction("flag", ("bounds",)),
            "clean": Prediction("pass"),
            "weak": Prediction("flag", ("snr",)),
            "variant": Prediction("flagged", ("satellites",)),
        }
        manifest = {"cases": cases, "variants": [variant], "candidates": [{"path": "candidate"}]}
        report = score_integrals(manifest, predictions, processing_version="test-v2")
        cohorts = report["cohorts"]
        self.assert_counts(cohorts["real"], labels=(2, 4, 1), statuses=(4, 1, 1, 1), flags=(1, 2, 1),
                           errors=(1, 2, 3), scores=(1/3, 1/4, 5/7, 4/6), assessed=(5, 4), groups=6)
        self.assert_counts(cohorts["synthetic"], labels=(1, 1, 0), statuses=(1, 1, 0, 0), flags=(0, 1, 0),
                           errors=(1, 0, 0), scores=(1, 1, 1, 1), assessed=(2, 2), groups=2)
        for name in ("real_processing_only", "variant"):
            self.assert_counts(cohorts[name], labels=(0, 1, 0), statuses=(1, 0, 0, 0), flags=(0, 1, 0),
                               errors=(1, 0, 0), scores=(1, 1, 1, 1), assessed=(1, 1), groups=1)
        self.assert_causes(cohorts["real"], [
            ("baseline", 2, 1, 0, 0, "test-v2"), ("bounds", 1, 1, 1, 0, "test-v2"),
            ("low_snr", 0, 1, 0, 0, None), ("neighbour", 0, 0, 0, 1, None),
            ("phase", 1, 1, 1, 0, "test-v2"), ("satellites", 1, 0, 0, 0, None),
        ])
        self.assert_causes(cohorts["synthetic"], [("low_snr", 1, 1, 1, 0, None)])
        self.assert_causes(cohorts["variant"], [("satellites", 1, 1, 1, 0, None)])
        self.assert_causes(cohorts["real_processing_only"], [("bounds", 1, 1, 1, 0, "test-v2")])
        self.assertEqual(report["groups"], {
            "hit": ["hit", "wrong", "control"], "unclear": ["unclear"], "missing": ["missing"],
            "false_flag": ["false_flag"], "pass": ["pass"], "unknown": ["unknown"],
            "clean": ["clean", "variant"], "weak": ["weak"],
        })
        self.assertEqual(len(report["rows"]), 11)
        self.assertEqual(report["rows"][7]["cohort"], "real_processing_only")
        self.assertEqual(report["candidates"], 1)
        self.assertIsNone(report["processing_control_note"])

    def test_zero_denominators_and_no_identified_controls(self):
        manifest = {"cases": [row("unknown", "indeterminate", "neighbour")], "variants": [], "candidates": []}
        report = score_integrals(manifest, {"unknown": Prediction("pass")}, processing_version="test-v2")
        self.assert_counts(report["cohorts"]["real"], labels=(0, 0, 1), statuses=(0, 1, 0, 0), flags=(0, 0, 0),
                           errors=(0, 0, 0), scores=(None, None, 1, None), assessed=(1, 0), groups=1)
        for name in ("synthetic", "variant", "real_processing_only"):
            self.assert_counts(report["cohorts"][name], labels=(0, 0, 0), statuses=(0, 0, 0, 0), flags=(0, 0, 0),
                               errors=(0, 0, 0), scores=(None, None, None, None), assessed=(0, 0), groups=0)
        self.assertEqual(report["processing_control_note"], "No real processing-only controls identified")

    def test_frozen_control_id_requires_the_same_source_and_bounds(self):
        control = next(case for case in load_cases()["cases"] if case["id"] == "bruker_1_15")
        for change in ({"source": "other/10"}, {"bounds": [1.14, 1.22]}):
            with self.subTest(change=change):
                manifest = {"cases": [{**control, **change}], "variants": [], "candidates": []}
                predictions = {control["id"]: Prediction("flag", ("bounds",))}
                report = score_integrals(manifest, predictions, processing_version="test-v2")
                self.assertEqual(report["rows"][0]["cohort"], "real")
                self.assertEqual(report["cohorts"]["real"]["true_positive"], 1)
                self.assertEqual(report["cohorts"]["real_processing_only"]["labelled"], 0)

    def test_frozen_fixture_counts(self):
        self.assertEqual(self.fixture["processing_version"], PROCESSING_VERSION)
        rows = {case["id"]: case for case in self.fixture["rows"]}
        for name, status, causes in (
            ("foreign_line", "flag", ["neighbour"]),
            ("foreign_tail", "flag", ["bounds", "neighbour"]),
            ("no_noise_window", "indeterminate", []),
        ):
            self.assertEqual((rows[name]["status"], rows[name]["flag_causes"]), (status, causes), name)
        cohorts = self.fixture["cohorts"]
        # FID02 is flagged twice on unknown truth. The absent FID05 is a miss.
        self.assert_counts(cohorts["real"], labels=(0, 1, 2), statuses=(2, 0, 0, 1), flags=(0, 0, 2),
                           errors=(0, 0, 1), scores=(None, 0, 2/3, 0), assessed=(2, 0), groups=2)
        # Six unreliable rows match. The low S/N row flags baseline only and
        # the local baseline row passes. Both indeterminate rows abstain.
        self.assert_counts(cohorts["synthetic"], labels=(3, 8, 2), statuses=(7, 4, 2, 0), flags=(0, 7, 0),
                           errors=(6, 1, 2), scores=(6/7, 6/8, 11/13, 11/11), assessed=(11, 11), groups=8)
        self.assert_counts(cohorts["variant"], labels=(0, 1, 0), statuses=(1, 0, 0, 0), flags=(0, 1, 0),
                           errors=(1, 0, 0), scores=(1, 1, 1, 1), assessed=(1, 1), groups=1)
        self.assert_counts(cohorts["real_processing_only"], labels=(0, 2, 0), statuses=(2, 0, 0, 0), flags=(0, 2, 0),
                           errors=(2, 0, 0), scores=(1, 1, 1, 1), assessed=(2, 2), groups=1)
        self.assert_causes(cohorts["real"], [("baseline", 0, 0, 0, 2, PROCESSING_VERSION),
                                            ("neighbour", 1, 0, 0, 0, None)])
        self.assert_causes(cohorts["synthetic"], [
            ("baseline", 1, 2, 0, 0, PROCESSING_VERSION), ("bounds", 3, 3, 3, 0, PROCESSING_VERSION),
            ("capture", 1, 3, 1, 0, None), ("low_snr", 1, 0, 0, 0, None), ("neighbour", 2, 2, 2, 0, None),
            ("phase", 1, 2, 1, 0, PROCESSING_VERSION), ("satellites", 1, 1, 1, 0, None),
        ])
        self.assert_causes(cohorts["variant"], [("satellites", 1, 1, 1, 0, None)])
        self.assert_causes(cohorts["real_processing_only"], [
            ("baseline", 1, 2, 1, 0, PROCESSING_VERSION), ("bounds", 2, 2, 2, 0, PROCESSING_VERSION),
            ("phase", 0, 2, 0, 0, PROCESSING_VERSION),
        ])
        self.assertEqual(self.fixture["candidates"], 3)
        self.assertEqual(len(self.fixture["rows"]), 19)
        groups = self.fixture["groups"]
        self.assertEqual(len(groups), 11)
        self.assertEqual(groups["bruker_1_15"], ["bruker_1_15", "bruker_7_64"])
        self.assertEqual(groups["fid02_dmso2"], ["fid02_dmso2", "fid02_ma"])
        self.assertEqual(groups["satellites_inside"], ["satellites_inside", "satellites_outside", "one_satellite_clipped"])

    def test_checker_refuses_mismatched_synthetic_source(self):
        rows = [case for case in load_cases()["cases"] if case["id"] in ("clean_singlet", "assigned_triplet")]
        prepared, _ = prepare_fixture(rows)
        singlet = prepared["clean_singlet"]
        result = check_region(singlet["spectrum"], singlet["bounds"], singlet["components"],
                              raw=prepared["assigned_triplet"]["raw"], noise_bounds_ppm=singlet["noise_bounds_ppm"])
        self.assertEqual(result.status, "cannot be assessed")
        self.assertEqual(result.snr.reason, "Raw source does not match the processed spectrum")

    def test_scores_are_independent_of_evaluator_input_annotations(self):
        truth = load_cases()
        rows = (*truth["cases"], *truth["variants"])
        row_ids = {row["id"] for row in rows}
        annotations_only = {"label", "causes", "evidence", "checker_status", "purity_eligible", "edge_step_ppm",
                            "edge_area_relative_changes", "half_block_relative_change", "area_cv"}

        def evaluate(inputs):
            self.assertEqual(len(inputs), len(rows) - 1)  # FID05 is absent.
            keys = {case["key"] for case in inputs}
            self.assertEqual(len(keys), len(inputs))
            self.assertFalse(keys & row_ids)
            for case in inputs:
                self.assertEqual(case.keys(), {"key", "raw", "spectrum", "bounds", "components", "noise_bounds_ppm"})
                if case["raw"].vendor == "synthetic":
                    self.assertRegex(case["raw"].source.name, r"^synthetic-[0-9a-f]{32}$")
                    self.assertEqual(case["spectrum"].source, case["raw"].source)
            return fixture_predictions(inputs)

        for annotations in ({}, dict.fromkeys(annotations_only, "scrambled")):
            with self.subTest(annotations=annotations):
                cases = [{**{k: v for k, v in case.items() if k not in annotations_only}, **annotations} for case in rows]
                prepared, version = prepare_fixture(cases)
                self.assertEqual(prepared["one_satellite_clipped"]["raw"].source,
                                 prepared["satellites_inside"]["raw"].source)
                predictions = evaluate_integrals(prepared, evaluate)
                self.assertEqual(score_integrals(truth, predictions, processing_version=version), self.fixture)


if __name__ == "__main__":
    unittest.main()
