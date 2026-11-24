"""Compare the six Nishizaki pilot FIDs with the workbook areas in docs/pilot-pairs.md.

From the repository, run python tools/compare_references.py DATA_ROOT.
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from dataclasses import replace
from math import isfinite
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
REFERENCES_FILE = PROJECT / "tests/fixtures/processing_references.json"
INTEGRAL_CASES_FILE = PROJECT / "tests/fixtures/integral_cases.json"
DEPOSIT = "334_Dataverse_20200705"
EXTRACTED = Path(f"x_{DEPOSIT}") / DEPOSIT


def load_references():
    return json.loads(REFERENCES_FILE.read_text(encoding="utf-8"))


def _pilot_paths(pilot):
    if len(pilot["runs"]) != 2:
        raise ValueError("Pilot requires two runs")
    paths = {}
    for run in pilot["runs"]:
        for cell in run["regions"]:
            if cell in paths:
                raise ValueError(f"Pilot cell {cell} appears in both runs")
            paths[cell] = run["path"]

    roles = set()
    for numerator, denominator in pilot.get("ratios", ()):
        sample, ec = paths[numerator], paths[denominator]
        if sample != ec:
            roles.add((sample, ec))
    if not roles:
        raise ValueError("Pilot requires a cross-run ratio")
    if len(roles) != 1:
        raise ValueError("Cross-run ratios disagree on sample and EC roles")
    return roles.pop()


def measure_regions(path, regions_ppm, references, *, processing_override=None, **reference):
    """Process all declared regions together and measure them on the corrected trace."""
    from spinivo import read
    from spinivo.measurement import measure_region
    from spinivo.processing import ProcessingSettings, process

    raw = read(path)
    noise = tuple(references["noise_windows_ppm"][raw.vendor])
    options = {**references["processing"], **(processing_override or {}), **reference}
    settings = ProcessingSettings(**options)
    spectrum = process(raw, settings, regions_ppm=tuple(map(tuple, regions_ppm.values())))
    return {name: measure_region(spectrum, name, tuple(bounds), noise)
            for name, bounds in regions_ppm.items()}


def report_pilot(root, pilot, references):
    regions, areas, failures, paths = {}, {}, {}, {}
    for run in pilot["runs"]:
        path = root / run["path"]
        regions.update(run["regions"])
        for cell in run["regions"]:
            paths[cell] = path
        try:
            measured = measure_regions(path, {cell: region["bounds_ppm"] for cell, region in run["regions"].items()},
                                       references)
        except (OSError, ValueError) as error:
            failures.update(dict.fromkeys(run["regions"], str(error)))
            continue
        for cell, region in measured.items():
            if region.area is not None and isfinite(region.area):
                areas[cell] = region.area
            else:
                failures[cell] = region.area_reason or "nonfinite area"

    for numerator, denominator in pilot["ratios"]:
        top, bottom = regions[numerator], regions[denominator]
        reference = top["area"] / bottom["area"]
        gaps = ["workbook integration bounds and intensity scaling undocumented"]
        if paths[numerator] != paths[denominator]:
            gaps.append("cross-run response and receiver scaling unverified")
        gaps.extend(record["evidence_gap"] for record in (pilot, top, bottom) if "evidence_gap" in record)
        gaps.extend(f"{cell}: {failures[cell]}" for cell in (numerator, denominator) if cell in failures)

        ours = residual = "unavailable"
        if numerator in areas and denominator in areas:
            if areas[denominator] == 0:
                gaps.append(f"{denominator}: zero denominator area")
            else:
                ratio = areas[numerator] / areas[denominator]
                ours, residual = f"{ratio:.8g}", f"{(ratio - reference) / reference:+.4%}"
                if min(areas[numerator], areas[denominator]) <= 0:
                    gaps.append("nonpositive signal area")
        print(f'{pilot["name"]} {top["signal"]}/{bottom["signal"]} ({numerator}/{denominator}) '
              f'reference={reference:.8g} ours={ours} relative_residual={residual} evidence_gaps={"; ".join(gaps)}')


def _quantitative_details(spectrum, result):
    def percent(value, fmt=".4%"):
        return "unavailable" if value is None or not isfinite(value) else format(value, fmt)

    causes = ",".join(dict.fromkeys(check for check, _ in result.flags)) or "none"
    recovery = result.baseline.details.recovery
    position = recovery.position_sensitivity
    position_pp = "unavailable" if position is None or not isfinite(position) else f"{100 * position:.4f}pp"
    fit = spectrum.smooth_baseline_fit
    half_block = next((region.half_block_relative_change for region in fit.regions
                       if region.bounds_ppm == result.bounds_ppm), None) if fit else None
    edges = " ".join(f"{trial.move}={percent(trial.relative_change, '+.4%')}" for trial in result.bounds.details)
    return (f"  causes={causes} threshold_version={result.thresholds.version} "
            f"midpoint_recovery={percent(recovery.value)} P={position_pp} "
            f"C={percent(recovery.placement_sensitivity)} half_block_change={percent(half_block, '+.4%')}\n"
            f"  edge_moves: {edges}")


def report_checks(root, references, *, quantitative=False):
    from spinivo import read
    from spinivo.integral_checks import PILOT
    from spinivo.integral_components import ComponentEvidence, check_region
    from spinivo.integral_report import region_report, render_text
    from spinivo.processing import QUANTITATIVE_SETTINGS, ProcessingSettings, process

    manifest = json.loads(INTEGRAL_CASES_FILE.read_text(encoding="utf-8"))
    defaults = QUANTITATIVE_SETTINGS if quantitative else ProcessingSettings(**references["rubric_processing"])
    runs = defaultdict(list)
    for case in manifest["cases"]:
        if case["kind"] != "real":
            continue
        source = Path(case["source"])
        path = root / source.relative_to(DEPOSIT) if source.parts[0] == DEPOSIT else PROJECT / source
        if path.exists():
            shift = case.get("reference_shift_ppm", defaults.reference_shift_ppm)
            runs[path.resolve(), shift].append(case)

    reports = {}
    for (path, shift), cases in runs.items():
        raw = read(path)
        reference_source = next((case["reference_source"] for case in cases if "reference_source" in case),
                                defaults.reference_source)
        settings = replace(defaults, reference_shift_ppm=shift, reference_source=reference_source)
        regions = tuple(tuple(case["bounds"]) for case in cases)
        spectrum = process(raw, settings, regions_ppm=regions)
        noise = tuple(references["noise_windows_ppm"][raw.vendor])
        for case, bounds in zip(cases, regions):
            result = check_region(spectrum, bounds, ComponentEvidence(), raw=raw, noise_bounds_ppm=noise,
                                  thresholds=PILOT)
            report = region_report(result, run=case["source"], signal=case["id"], thresholds=PILOT)
            reports[case["id"]] = render_text(report)
            if quantitative:
                reports[case["id"]] += "\n" + _quantitative_details(spectrum, result)

    for case in manifest["cases"]:
        if case["id"] in reports:
            print(reports[case["id"]])
            print(f'label={case["label"]}')


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("data_root", type=Path, help="local data cache or extracted Nishizaki deposit directory")
    parser.add_argument("--checks", action="store_true", help="report checks for available real integral cases")
    parser.add_argument("--quantitative", action="store_true", help="use the quantitative profile with --checks")
    args = parser.parse_args()
    if args.quantitative and not args.checks:
        parser.error("--quantitative requires --checks")
    root = args.data_root
    if not root.is_dir():
        parser.error(f"data root is not a directory: {root}")
    for relative in (Path("nishizaki_QFY6L9") / EXTRACTED, EXTRACTED, Path(DEPOSIT)):
        if (root / relative).is_dir():
            root /= relative
            break
    references = load_references()
    if args.checks:
        report_checks(root, references, quantitative=args.quantitative)
    else:
        from tools.build_case_manifest import _nishizaki_key
        from tools.check_case_policy import preflight_purity

        for pilot in references["pilots"]:
            try:
                sample, ec = map(_nishizaki_key, _pilot_paths(pilot))
                preflight_purity("external_calibration", sample, ec)
            except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
                detail = " ".join(str(error).splitlines())
                parser.exit(1, f"Pilot {pilot['name']} refused: {detail}\n")
        for pilot in references["pilots"]:
            report_pilot(root, pilot, references)


if __name__ == "__main__":
    sys.path.insert(0, str(PROJECT))
    sys.path.insert(0, str(PROJECT / "src"))
    main()
