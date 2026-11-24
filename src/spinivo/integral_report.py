"""Render integral region checks as versioned JSON and plain text."""

from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from math import isfinite
from numbers import Integral, Real

from .integral_checks import IntegralThresholds
from .integral_components import ComponentObservation, RegionResult


def _json_value(value):
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        number = float(value)
        return number if isfinite(number) else None
    if is_dataclass(value):
        value = asdict(value)
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    raise TypeError(f"Unsupported report value: {type(value).__name__}")


def _observation_report(observation):
    if isinstance(observation, ComponentObservation):
        details = {"bands": observation.bands} if observation.check == "satellites" else None
    else:
        details = observation.details
    return {
        "check": observation.check, "status": observation.status,
        "value": observation.value, "unit": observation.unit, "limit": observation.limit,
        "reason": observation.reason, "threshold_version": observation.threshold_version, "details": details,
    }


def _threshold_report(result, thresholds):
    """Verify the supplied profile and recorded observations against the profile used at evaluation."""
    if thresholds not in (result.thresholds, result.requested_thresholds):
        raise ValueError("Threshold profile does not match the region observations")
    profile = asdict(result.thresholds)
    limits = {"area": None, "snr": profile["min_snr"]}
    if any(observation.threshold_version != profile["version"]
           or observation.limit != limits.get(observation.check, profile["max_relative_change"])
           for observation in result.observations):
        raise ValueError("Threshold profile does not match the region observations")
    return profile


def region_report(result: RegionResult, *, run: str, signal: str, thresholds: IntegralThresholds) -> dict:
    """Build a JSON-compatible report with the profile used to check the region."""
    evidence = result.evidence
    return _json_value({
        "schema": "spinivo.integral_region/1",
        "run": run,
        "signal": signal,
        "bounds_ppm": result.bounds_ppm,
        "nominal_area": {"value": result.nominal_area, "unit": result.area.unit},
        "status": result.status,
        "purity_eligible": result.status == "pass",
        "flags": [{"check": check, "reason": reason} for check, reason in result.flags],
        "gaps": [{"check": check, "reason": reason} for check, reason in result.gaps],
        "thresholds": _threshold_report(result, thresholds),
        "evidence": {"target_source": evidence.target_source, "foreign_source": evidence.foreign_source,
                     "satellite_source": evidence.satellite_source, "satellite_convention": evidence.satellite_convention},
        "observations": [_observation_report(observation) for observation in result.observations],
    })


def render_json(report: dict) -> str:
    """Render a region report with deterministic keys and strict JSON numbers."""
    return json.dumps(report, sort_keys=True, allow_nan=False)


def _bounds_text(bounds):
    return "-".join(f"{value:.12g}" if value is not None else "unavailable" for value in bounds)


def render_text(report: dict) -> str:
    """Render the region outcome and its observations."""
    bounds = _bounds_text(report["bounds_ppm"])
    lines = [f'{report["signal"]} {report["run"]} {bounds} ppm: {report["status"].upper()}']
    if report["status"] == "pass":
        area = report["nominal_area"]
        value = f'{area["value"]:.8g}' if area["value"] is not None else "unavailable"
        lines.append(f'  Area: {value} {area["unit"]}')

    for flag in report["flags"]:
        lines.append(f'  {flag["check"]}: {flag["reason"]}')
        if flag["check"] != "bounds":
            continue
        bounds_observation = next(obs for obs in report["observations"] if obs["check"] == "bounds")
        for trial in bounds_observation["details"]:
            change = trial["relative_change"]
            change_text = f"{change:+.4%}" if change is not None else f'unavailable ({trial["reason"]})'
            lines.append(f'    {trial["move"]} {_bounds_text(trial["bounds_ppm"])} ppm: {change_text}')

    if report["gaps"]:
        lines.append("  Evidence gaps")
        for gap in report["gaps"]:
            lines.append(f'    {gap["check"]}: {gap["reason"]}')

    baseline = next(obs for obs in report["observations"] if obs["check"] == "baseline")
    flips = baseline["details"]["recovery"].get("anchor_flips", ())
    if flips:
        lines.append("  Spike anchor changes without the lock (observations only)")
        for trial in flips:
            for key, label in (("added_ppm", "would add"), ("removed_ppm", "would remove")):
                if trial[key]:
                    centers = ", ".join(f"{ppm:.6g}" for ppm in trial[key])
                    lines.append(f'    {trial["position_percent"]}%: {label} {centers} ppm')
    return "\n".join(lines)
