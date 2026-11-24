"""Score frozen integral labels against independent regional checker results."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from uuid import uuid4


_RECIPE = {"target": "singlet", "foreign": None, "phase_deg": 0, "anchor_scale": 1,
           "baseline_area": 0, "snr": 3000}
_EVALUATOR_FIELDS = ("raw", "spectrum", "bounds", "components", "noise_bounds_ppm")
_COHORTS = ("real", "synthetic", "variant", "real_processing_only")
# The frozen pilot rubric identifies these two windows as processing controls.
_RUBRIC_CONTROL_SOURCE = "tests/fixtures/bruker/zenodo17083860_833-1/10"
_RUBRIC_CONTROLS = {"bruker_1_15": (1.15, 1.22), "bruker_7_64": (7.64, 7.69)}
_PROCESSING_CAUSES = {"baseline", "bounds", "phase"}
_LABELS = ("reliable", "unreliable", "indeterminate")
_STATUSES = {"flag": "flag", "flagged": "flag", "pass": "pass", "indeterminate": "indeterminate",
             "cannot be assessed": "indeterminate", "not assessed": "not assessed"}


@dataclass(frozen=True)
class Prediction:
    """A regional outcome and the check names that supported its flags."""

    status: str
    flags: tuple[str, ...] = ()

    @classmethod
    def from_region(cls, result):
        return cls(result.status, tuple(check for check, _ in result.flags))


def evaluate_integrals(prepared, evaluate):
    """Evaluate prepared traces and declared checker inputs keyed by row ID.

    The evaluator sees randomized keys and returns predictions keyed by them.
    """
    by_key = {uuid4().hex: row_id for row_id in prepared}
    inputs = tuple({"key": key, **{name: prepared[by_key[key]][name] for name in _EVALUATOR_FIELDS}}
                   for key in sorted(by_key))
    return {by_key[key]: prediction for key, prediction in evaluate(inputs).items()}


def _source_group(row, by_id):
    source = row["kind"], row["source"]
    if row["kind"] == "synthetic":
        # The generator file is shared by distinct FIDs. A repeated recipe at
        # the nominal seed and a changed window on its parent share one trace.
        trace = by_id[row["parent"]] if "parent" in row else row
        return (*source, *(trace.get(key, default) for key, default in _RECIPE.items()))
    return source


def _scored_row(row, cohort, group, prediction):
    if row["label"] not in _LABELS:
        raise ValueError(f"Unknown label for {row['id']}: {row['label']}")
    prediction = prediction or Prediction("not assessed")
    if prediction.status not in _STATUSES:
        raise ValueError(f"Unknown prediction for {row['id']}: {prediction.status}")
    status = _STATUSES[prediction.status]
    flagged = sorted({"low_snr" if check == "snr" else check for check in prediction.flags})
    if flagged and status != "flag":
        raise ValueError(f"Flags require a flagged region: {row['id']}")
    matched = sorted(set(row["causes"]).intersection(flagged)) if row["label"] == "unreliable" else []
    return {**row, "cohort": cohort, "group": group, "status": status,
            "flag_causes": flagged, "matched_causes": matched}


def _fraction(count, total):
    return count / total if total else None


def _summarize(rows, processing_version):
    labels = Counter(row["label"] for row in rows)
    statuses = Counter(row["status"] for row in rows)
    determinate = [row for row in rows if row["label"] != "indeterminate"]
    flags = Counter(row["label"] for row in rows if row["status"] == "flag")
    assessed = statuses["flag"] + statuses["pass"]
    determinate_assessed = sum(row["status"] in ("flag", "pass") for row in determinate)
    determinate_flags = flags["reliable"] + flags["unreliable"]
    tp = sum(bool(row["matched_causes"]) for row in determinate)

    expected = Counter(cause for row in determinate if row["label"] == "unreliable" for cause in set(row["causes"]))
    detected = Counter(cause for row in determinate for cause in row["flag_causes"])
    matched = Counter(cause for row in determinate for cause in row["matched_causes"])
    indeterminate_flags = Counter(cause for row in rows if row["label"] == "indeterminate"
                                  for cause in row["flag_causes"])
    causes = [
        {"cause": cause, "expected": expected[cause], "flagged": detected[cause], "matched": matched[cause],
         "flagged_indeterminate": indeterminate_flags[cause],
         "processing_version": processing_version if cause in _PROCESSING_CAUSES else None}
        for cause in sorted(expected.keys() | detected.keys() | indeterminate_flags.keys())
    ]
    return {
        "labelled": len(rows), "determinate": len(determinate), "groups": len({row["group"] for row in rows}),
        "labels": {label: labels[label] for label in _LABELS},
        "statuses": {status: statuses[status] for status in ("flag", "pass", "indeterminate", "not assessed")},
        "region_flags": {label: flags[label] for label in _LABELS},
        "true_positive": tp, "false_positive": determinate_flags - tp, "false_negative": labels["unreliable"] - tp,
        "precision": _fraction(tp, determinate_flags), "recall": _fraction(tp, labels["unreliable"]),
        "assessed": assessed, "determinate_assessed": determinate_assessed,
        "coverage": _fraction(assessed, len(rows)), "determinate_coverage": _fraction(determinate_assessed, len(determinate)),
        "causes": causes,
    }


def score_integrals(manifest, predictions, *, processing_version):
    """Summarize row-level scores and identify correlated rows by source.

    A flagged unreliable row needs at least one matching cause to be a hit.
    A wrong-cause flag on a known unreliable row counts as both a false
    positive and a false negative. Missing predictions count as abstentions.
    Cause counts use determinate labels, with flags on indeterminate labels
    recorded separately. Source groups describe correlation in the row counts.
    """
    cases, variants = manifest["cases"], manifest["variants"]
    by_id = {row["id"]: row for row in (*cases, *variants)}
    if len(by_id) != len(cases) + len(variants):
        raise ValueError("Integral row IDs must be unique")
    cohort_rows = {cohort: [] for cohort in _COHORTS}
    group_ids, groups, rows = {}, defaultdict(list), []
    for cohort, row_set in ((None, cases), ("variant", variants)):
        for row in row_set:
            name = cohort or row["kind"]
            if name == "real" and (row.get("processing_only") or (
                row["source"] == _RUBRIC_CONTROL_SOURCE and tuple(row["bounds"]) == _RUBRIC_CONTROLS.get(row["id"])
            )):
                name = "real_processing_only"
            key = _source_group(row, by_id)
            group = group_ids.setdefault(key, row["id"])
            groups[group].append(row["id"])
            scored = _scored_row(row, name, group, predictions.get(row["id"]))
            cohort_rows[name].append(scored)
            rows.append(scored)

    return {
        "processing_version": processing_version,
        "cohorts": {name: _summarize(items, processing_version) for name, items in cohort_rows.items()},
        "rows": rows, "groups": dict(groups),
        "candidates": len(manifest["candidates"]),
        "processing_control_note": (None if cohort_rows["real_processing_only"]
                                    else "No real processing-only controls identified"),
    }


def score_fixture(path, evaluate, *, prepare):
    """Read the frozen fixture manifest and evaluate its labelled rows."""
    manifest = json.loads(Path(path).read_text(encoding="utf-8"))
    prepared, processing_version = prepare((*manifest["cases"], *manifest["variants"]))
    predictions = evaluate_integrals(prepared, evaluate)
    return score_integrals(manifest, predictions, processing_version=processing_version)
