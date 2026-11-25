"""Check the frozen baseline cohorts and enforce the purity pilot allowlist."""

import argparse
import csv
import json
import sys
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.check_case_manifest import INVENTORY, MANIFEST, PROJECT, validate

WINDOWS = PROJECT / "data/baseline_windows.json"
WIDENED = "fid05_widened_dmso2"


def check_cohorts(manifest_path=MANIFEST, windows_path=WINDOWS):
    """Return the recomputed cohort records after checking the snapshot and membership."""
    from tools.build_case_manifest import _window_cases

    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    snapshot = json.loads(Path(windows_path).read_text(encoding="utf-8"))
    widened = [case for case in manifest["cases"] if case["id"] == WIDENED]
    if len(widened) != 1 or widened[0]["kind"] != "overlap" or widened[0]["cohorts"] != []:
        raise ValueError(f"{WIDENED} must be an overlap case outside the baseline cohorts")

    windows = []
    for case in manifest["cases"]:
        if case["kind"] == "window":
            windows.append(case)
        elif case.get("cohorts"):
            raise ValueError(f"Non-window case {case['id']} declares baseline cohort membership")

    expected, cohorts = _window_cases(snapshot)
    computed, recomputed = _window_cases(windows)
    if sorted(case["id"] for case in computed) != sorted(case["id"] for case in expected):
        raise ValueError("Baseline windows differ from the frozen snapshot")
    for case, calculated in zip(windows, computed):
        if case["id"] != calculated["id"] or case["cohorts"] != calculated["cohorts"]:
            raise ValueError(f"Incorrect baseline window ID or cohort membership for {case['id']}")
    if recomputed != cohorts or manifest["cohorts"] != cohorts:
        raise ValueError("Baseline cohort records differ from the frozen snapshot")
    return recomputed


def _check_pilot_pairs(manifest):
    from tools.build_case_manifest import PILOT_PATHS, _nishizaki_key

    pins = {f"pilot-{figure}": tuple(map(_nishizaki_key, paths)) for figure, paths in PILOT_PATHS.items()}
    pairs = {pair["id"]: (tuple(pair["sample"]["raw_key"]), tuple(pair["ec"]["raw_key"]))
             for pair in manifest["pilot_pairs"]}
    if pairs != pins or len(manifest["pilot_pairs"]) != len(pins):
        raise ValueError("Manifest pilot pairs differ from the pinned pilot pairs")
    return set(pairs.values())


def preflight_purity(method, sample_key, ec_key=None):
    """Require eligible inventory keys before reading a real purity pilot."""
    components = validate(MANIFEST, INVENTORY)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    pairs = _check_pilot_pairs(manifest)
    allowed = {key for pair in pairs for key in pair}
    if method not in ("internal_standard", "external_calibration"):
        raise ValueError(f"Unknown purity method: {method!r}")
    if method == "internal_standard" and ec_key is not None:
        raise ValueError("Internal standard purity takes a sample key only")
    keys = []
    for role, key in (("sample", sample_key), ("EC", ec_key)):
        if role == "EC" and method == "internal_standard":
            continue
        if (not isinstance(key, (tuple, list)) or not all(isinstance(part, str) for part in key)
                or tuple(key) not in allowed):
            raise ValueError(f"{role} key is not a pinned pilot key: {key!r}")
        keys.append(tuple(key))
    if method == "external_calibration" and tuple(keys) not in pairs:
        raise ValueError("External calibration requires a pinned pilot pair in sample, EC order")

    runs = {run["id"] for run in manifest["runs"] if tuple(run["raw_key"]) in keys}
    holdouts = {("case", case["id"]) for case in manifest["cases"] if case["split"] == "final_holdout"}
    exposures = {("run", run_id) for exposure in manifest["exposures"] for run_id in exposure.get("run_ids") or ()}
    for component in components:
        if not any(kind == "run" and identity in runs for kind, identity in component):
            continue
        if holdouts.intersection(component):
            raise ValueError("Purity pilot belongs to a final_holdout component")
        if exposures.intersection(component):
            raise ValueError("Purity pilot belongs to an exposed component")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--inventory", type=Path, default=INVENTORY)
    parser.add_argument("--windows", type=Path, default=WINDOWS)
    args = parser.parse_args()
    try:
        validate(args.manifest, args.inventory)
        _check_pilot_pairs(json.loads(args.manifest.read_text(encoding="utf-8")))
        cohorts = check_cohorts(args.manifest, args.windows)
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
        detail = " ".join(str(error).splitlines())
        print(f"case policy: INVALID ({detail})", file=sys.stderr)
        return 1
    counts = ", ".join(f"{name}: {len(cohort['members'])}" for name, cohort in cohorts.items())
    print(f"case policy: OK ({counts})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
