"""Audit one-setting changes to the pinned Nishizaki pilot area ratios."""

import argparse
import csv
import sys
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from math import isfinite
from pathlib import Path

if __name__ == "__main__":
    project = Path(__file__).resolve().parents[1]
    sys.path[:0] = [str(project), str(project / "src")]

from spinivo import read
from spinivo.assessment_report import baseline_shape_text
from spinivo.filters import resolve_delay
from spinivo.internal_standard import _finite
from spinivo.pipeline import BaselineShapeRegion, baseline_shape_diagnostic
from spinivo.processing import ProcessingSettings, baseline_shape_areas, integrate, process, uses_quantitative_settings
from tools.build_case_manifest import _nishizaki_key
from tools.check_case_policy import preflight_purity
from tools.compare_references import DEPOSIT, EXTRACTED, _pilot_paths, load_references


@dataclass(frozen=True)
class Settings:
    zero_fill_factor: int
    phase_lb_hz: float | None
    scan_normalisation: bool = False
    gain_normalisation: bool = False
    undo_reader_scale: bool = False
    remove_filter_delay: bool = True
    region_placement: str = "pilot"
    pulse_direction: str = "none"


def stages(processing):
    baseline = Settings(processing.zero_fill_factor, processing.phase_lb_hz)
    yield "baseline", "as in compare_references", baseline
    # Keep a phase copy in both cases so the baseline uses the same fit policy.
    for stage, setting, value, detail in (
        ("scan count", "scan_normalisation", True, "divide by NS if stored values are sums"),
        ("receiver gain", "gain_normalisation", True, "receiver response correction unknown"),
        ("reader scaling", "undo_reader_scale", True, "divide by documented reader multiplier"),
        ("filter", "remove_filter_delay", False, "bypass JEOL group delay correction"),
        ("axis sampling", "zero_fill_factor", 1, "zero fill factor 1, same header ppm reference"),
        ("phase", "phase_lb_hz", processing.lb_hz, f"fit on LB{processing.lb_hz:g} copy"),
        ("region placement", "region_placement", "sample", "move sample upper edge outward"),
        ("region placement", "region_placement", "EC", "move EC upper edge outward"),
        ("region placement", "region_placement", "workbook", "workbook bounds unknown"),
        ("pulse direction", "pulse_direction", "sample/EC", "nominal width sample/EC, arithmetic only"),
        ("pulse direction", "pulse_direction", "EC/sample", "nominal width EC/sample, arithmetic only"),
    ):
        yield stage, detail, replace(baseline, **{setting: value})


def quotient(top, bottom):
    if top is None or bottom is None or not isfinite(top) or not isfinite(bottom) or not bottom:
        return None
    result = top / bottom
    return result if isfinite(result) else None


def positive(value):
    return value if value is not None and isfinite(value) and value > 0 else None


def shown(value, fmt=".8g"):
    return "unknown" if value is None or not isfinite(value) else format(value, fmt)


def trial_areas(settings, area, raw, top, bottom):
    if settings.gain_normalisation or settings.region_placement == "workbook":
        return None, None, None
    areas = []
    for role, cell in (("sample", top), ("EC", bottom)):
        value = area(cell, settings, settings.region_placement == role)
        if settings.scan_normalisation:
            value = quotient(value, positive(raw[cell].scans))
        if settings.undo_reader_scale:
            value = quotient(value, positive(raw[cell].intensity_scale))
        areas.append(value)

    ratio = quotient(*areas)
    if settings.pulse_direction != "none":
        first, second = (top, bottom) if settings.pulse_direction == "sample/EC" else (bottom, top)
        factor = quotient(positive(raw[first].pulse_width_us), positive(raw[second].pulse_width_us))
        ratio = ratio * factor if ratio is not None and factor is not None else None
    return *areas, ratio


def report_pilot(root, pilot, references, paths):
    processing = ProcessingSettings(**references["processing"])
    runs = {run["path"]: run["regions"] for run in pilot["runs"]}
    raw_runs = {path: read(root / path) for path in paths}
    cells = {cell: (path, region) for path, regions in runs.items() for cell, region in regions.items()}
    raw = {cell: raw_runs[path] for cell, (path, _) in cells.items()}
    bounds = {}
    for cell, (_, region) in cells.items():
        limits = region.get("bounds_ppm")
        if isinstance(limits, (list, tuple)) and len(limits) == 2 and all(map(_finite, limits)):
            bounds[cell] = limits

    @lru_cache(maxsize=None)
    def spectrum(path, remove_filter_delay, zero_fill_factor, phase_lb_hz):
        fid = raw_runs[path]
        if not remove_filter_delay:
            fid = replace(fid, digital_filter={"digital_filter": False})
        options = replace(processing, zero_fill_factor=zero_fill_factor, phase_lb_hz=phase_lb_hz)
        regions = tuple(tuple(bounds[cell]) for cell in runs[path] if cell in bounds)
        try:
            return process(fid, options, regions_ppm=regions), None
        except ValueError as error:
            return None, str(error)

    errors = {}

    @lru_cache(maxsize=None)
    def area(cell, settings, move_edge=False):
        if cell not in bounds:
            errors[cell, settings] = "pilot bounds unknown"
            return None
        path = cells[cell][0]
        trace, error = spectrum(path, settings.remove_filter_delay, settings.zero_fill_factor, settings.phase_lb_hz)
        if error is not None:
            errors[cell, settings] = error
            return None
        try:
            low, high = bounds[cell]
            if move_edge:
                # Same step as the integral bounds check, on the fixed baseline trace.
                high += max(trace.spectral_width_hz / trace.acquired_points / trace.observe_frequency_mhz,
                            0.05 * (high - low))
            return integrate(trace, (low, high)).area
        except ValueError as error:
            errors[cell, settings] = str(error)
            return None

    variants = tuple(stages(processing))
    baseline = variants[0][2]
    quantitative = uses_quantitative_settings(processing)

    def shape_region(cell, role):
        path = cells[cell][0]
        span = tuple(bounds[cell]) if cell in bounds else None
        areas = None
        if quantitative and span is not None:
            trace, _ = spectrum(path, baseline.remove_filter_delay, baseline.zero_fill_factor, baseline.phase_lb_hz)
            if trace is not None:
                areas = baseline_shape_areas(trace, span)
        return BaselineShapeRegion(path, role, cell, span, areas)

    purpose = {"FID05": "solvent-transfer challenge", "FID06": "calculation reproduction"}.get(pilot["name"])
    print(pilot["name"] + (f" ({purpose})" if purpose else ""))
    for role, path in zip(("sample", "EC"), paths):
        fid = raw_runs[path]
        print(f"  {role} NS={shown(fid.scans)} nominal_RG={shown(fid.receiver_gain)} "
              f"reader_multiplier={shown(fid.intensity_scale)} pulse_us={shown(fid.pulse_width_us)} "
              f"delay_points={shown(resolve_delay(fid).delay_points)} carrier_ppm={shown(fid.carrier_ppm)}")
    print(f"  baseline integrates LB{processing.lb_hz:g} with phase fit on LB{shown(processing.phase_lb_hz)} "
          f"and zero fill factor {processing.zero_fill_factor}")
    print("  accumulation convention, workbook intensity scale, workbook ppm shift and physical pulse transfer unknown")
    if "evidence_gap" in pilot:
        print(f"  {pilot['evidence_gap']}")

    controls = [(top, bottom) for top, bottom in pilot["ratios"]
                if cells[top][0] == cells[bottom][0] and cells[top][1]["signal"] == "MA"]
    print("  Per-run MA controls (baseline)" if controls else "  Per-run MA controls absent")
    for top, bottom in controls:
        role = "sample" if cells[top][0] == paths[0] else "EC"
        measured = quotient(area(top, baseline), area(bottom, baseline))
        reference = quotient(cells[top][1]["area"], cells[bottom][1]["area"])
        print(f"    {role} {top}/{bottom} ratio={shown(measured)} workbook_area_ratio={shown(reference)}")

    for top, bottom in pilot["ratios"]:
        if cells[top][0] == cells[bottom][0]:
            continue
        reference = quotient(cells[top][1]["area"], cells[bottom][1]["area"])
        base_ratio = trial_areas(baseline, area, raw, top, bottom)[2]
        relative = quotient(base_ratio, reference)
        print(f"  Cross-tube {top}/{bottom} workbook_area_ratio={shown(reference)} "
              f"baseline_vs_workbook={shown(None if relative is None else relative - 1, '+.4%')}")
        print(f'    pilot bounds (ppm) sample={bounds.get(top, "unknown")} EC={bounds.get(bottom, "unknown")}')
        for cell in (top, bottom):
            if "evidence_gap" in cells[cell][1]:
                print(f'    {cells[cell][1]["evidence_gap"]}')
        print("    stage | setting changed | sample area | EC area | ratio | change vs baseline")
        for stage, detail, settings in variants:
            sample, ec, ratio = trial_areas(settings, area, raw, top, bottom)
            reasons = [f"{cell}: {errors[cell, settings]}" for cell in (top, bottom) if (cell, settings) in errors]
            reasons.extend(f"{cell}: nonpositive signal area" for cell, value in ((top, sample), (bottom, ec))
                           if value is not None and value <= 0)
            if reasons:
                detail += " (" + ". ".join(reasons) + ")"
            relative = quotient(ratio, base_ratio)
            change = None if relative is None else relative - 1
            print(f"    {stage} | {detail} | {shown(sample)} | {shown(ec)} | {shown(ratio)} | {shown(change, '+.4%')}")

    for top, bottom in pilot["ratios"]:
        roles = ("numerator", "denominator")
        shape = baseline_shape_diagnostic(
            tuple(shape_region(cell, role) for cell, role in zip((top, bottom), roles)),
            roles, same_tube=cells[top][0] == cells[bottom][0], quantitative=quantitative,
        )
        print(f"  {top}/{bottom} {baseline_shape_text(asdict(shape))}")

    # A sensitivity or closer workbook match alone cannot identify the cause.
    # The committed pilot evidence does not establish any such explanation.
    print("  unresolved; EC validation gate stays closed")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path, help="local data cache or extracted Nishizaki deposit directory")
    root = parser.parse_args(argv).data_root
    if not root.is_dir():
        parser.error(f"data root is not a directory: {root}")
    for relative in (Path("nishizaki_QFY6L9") / EXTRACTED, EXTRACTED, Path(DEPOSIT)):
        if (root / relative).is_dir():
            root /= relative
            break

    references = load_references()
    pairs = []
    for pilot in references["pilots"]:
        try:
            paths = _pilot_paths(pilot)
            preflight_purity("external_calibration", *map(_nishizaki_key, paths))
        except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
            detail = " ".join(str(error).splitlines())
            parser.exit(1, f"Pilot {pilot['name']} refused: {detail}\n")
        pairs.append((pilot, paths))

    status = 0
    for pilot, paths in pairs:
        try:
            report_pilot(root, pilot, references, paths)
        except (OSError, ValueError) as error:
            print(f"Pilot {pilot['name']} unavailable ({error})", file=sys.stderr)
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
