"""Evaluate predeclared qNMR preparations from per-run results."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from math import comb, hypot, isfinite, sqrt, ulp
from statistics import mean, median, stdev
from sys import float_info

import numpy as np

from .external_calibration import _missing


ESTIMATOR = (
    "Within each method, vendor and D1 stratum, average the passing repeats with a "
    "measurement and reference for each eligible preparation. Report sample repeat SD "
    "separately. Weight preparations with results equally for mean signed error, "
    "median absolute error and fraction within tolerance. Count results over all "
    "predeclared eligible preparations and count independent units by dependency group."
)
# Method targets in percentage points (docs/decisions.md, decision 13).
_TARGETS_PP = {"internal_standard": 0.50, "external_calibration": 0.75}
_STATUSES = ("pass", "fail", "cannot be assessed")


@dataclass(frozen=True)
class Preparation:
    preparation_id: str
    method: str
    vendor: str
    d1_stratum: str


@dataclass(frozen=True)
class AccuracyRow:
    # Qualify the run ID by method if the same acquisition supports both methods.
    run_id: str
    preparation: Preparation
    status: str
    dependency_group: str
    measured_purity_pct: float | None = None
    reference_purity_pct: float | None = None
    measurement_u_pp: float | None = None
    reference_u_pp: float | None = None
    reference_coverage_basis: str | None = None


def _fraction(count, total):
    return count / total if total else None


def _nonnegative(value):
    return value is not None and isfinite(value) and value >= 0


def _standard_u(repeats, covariances):
    if not repeats:
        return None, "No eligible results"
    blocks = defaultdict(list)
    for rows in repeats:
        weight = 1 / (len(repeats) * len(rows))
        for row in rows:
            if not _nonnegative(row.measurement_u_pp) or not _nonnegative(row.reference_u_pp):
                return None, f"Missing or invalid standard uncertainty for {row.run_id}"
            basis = row.reference_coverage_basis
            if not isinstance(basis, str) or _missing(basis):
                return None, f"Unknown reference coverage basis for {row.run_id}"
            u = hypot(row.measurement_u_pp, row.reference_u_pp)
            variance = u * u
            if not isfinite(variance):
                return None, f"Nonfinite signed-error variance for {row.run_id}"
            blocks[row.dependency_group].append((row.run_id, weight, variance))

    variance = 0.0
    for block in blocks.values():
        weights = np.array([weight for _, weight, _ in block])
        matrix = np.diag([var for _, _, var in block])
        for i, (left, _, left_var) in enumerate(block):
            for j, (right, _, right_var) in enumerate(block[:i]):
                covariance = covariances.get(frozenset((left, right)))
                if covariance is None or not isfinite(covariance):
                    return None, f"Missing or invalid signed-error covariance for {left} and {right}"
                product = left_var * right_var
                # Take roots separately when multiplying the variances loses range.
                bound = (sqrt(product) if float_info.min <= product and isfinite(product)
                         else sqrt(left_var) * sqrt(right_var))
                # A perfectly correlated pair (covariance u1*u2) may exceed the rounded root by a few ULPs.
                if abs(covariance) > (bound + 4 * ulp(bound) if bound else 0.0):
                    return None, f"Signed-error covariance exceeds the variance bound for {left} and {right}"
                matrix[i, j] = matrix[j, i] = covariance
        scale = np.max(np.abs(matrix))
        if scale and np.linalg.eigvalsh(matrix / scale)[0] < -1e-12:
            return None, "Signed-error covariance matrix is not positive semidefinite"
        # A tiny negative quadratic form can result from floating-point rounding.
        variance += max(0.0, float(weights @ matrix @ weights))
    return sqrt(variance), None


def _preparation_result(preparation, rows, exclusion, k, covariances):
    counts = Counter(row.status for row in rows)
    usable = [] if exclusion is not None else [
        row for row in rows if row.status == "pass"
        and row.measured_purity_pct is not None and row.reference_purity_pct is not None
    ]
    references = {row.reference_purity_pct for row in usable}
    if len(references) > 1:
        raise ValueError(f"Conflicting reference purities for {preparation.preparation_id}")
    values = [row.measured_purity_pct for row in usable]
    measured = mean(values) if values else None
    reference = next(iter(references), None)
    error = measured - reference if usable else None
    standard_u, reason = _standard_u([usable] if usable else [], covariances)
    comparison_u = k * standard_u if k is not None and standard_u is not None else None
    return usable, {
        **asdict(preparation), "exclusion": exclusion,
        "dependency_groups": sorted({row.dependency_group for row in rows}),
        "statuses": {status: counts[status] for status in _STATUSES},
        "eligible_repeats": len(usable), "measured_purity_pct": measured,
        "reference_purity_pct": reference,
        "error_pp": error, "repeat_sd_pp": stdev(values) if len(values) > 1 else None,
        "error_standard_u_pp": standard_u, "error_u_reason": reason,
        "standard_u_ratio": abs(error) / comparison_u if comparison_u else None,
    }


def _summarize(results, repeats, tolerance, covariances):
    eligible = [result for result in results if result["exclusion"] is None]
    scored = [result for result in eligible if result["error_pp"] is not None]
    errors = [abs(result["error_pp"]) for result in scored]
    ratios = [result["standard_u_ratio"] for result in scored if result["standard_u_ratio"] is not None]
    standard_u, reason = _standard_u(repeats, covariances)
    return {
        "tolerance_pp": tolerance,
        "eligible_preparations": len(eligible), "result_preparations": len(scored),
        "eligible_dependency_groups": len({group for result in eligible for group in result["dependency_groups"]}),
        "result_dependency_groups": len({row.dependency_group for rows in repeats for row in rows}),
        "eligible_result_fraction": _fraction(len(scored), len(eligible)),
        "bias_pp": mean(result["error_pp"] for result in scored) if scored else None,
        "median_absolute_error_pp": median(errors) if errors else None,
        "agreement_fraction": _fraction(sum(error <= tolerance for error in errors), len(errors)),
        "standard_u_comparisons": len(ratios),
        "standard_u_agreement_fraction": _fraction(sum(ratio <= 1 for ratio in ratios), len(ratios)),
        "bias_standard_u_pp": standard_u, "bias_u_reason": reason,
    }


def _roster_runs(rows, preparations, exclusions, k, covariances):
    if k is not None and (not isfinite(k) or k <= 0):
        raise ValueError("k must be finite and positive")
    roster = {}
    for preparation in preparations:
        method, vendor, d1 = preparation.method, preparation.vendor, preparation.d1_stratum
        key = method, vendor, d1, preparation.preparation_id
        if method not in _TARGETS_PP or not all(part.strip() for part in key):
            raise ValueError(f"Invalid preparation declaration: {preparation}")
        if key in roster:
            raise ValueError(f"Duplicate preparation declaration: {key}")
        roster[key] = preparation
    repeats = {preparation: [] for preparation in roster.values()}
    if any(preparation not in repeats or not isinstance(reason, str) or _missing(reason)
           for preparation, reason in exclusions.items()):
        raise ValueError("Exclusions require a declared preparation and a reason")

    by_run = {}
    for row in rows:
        if (row.preparation not in repeats or row.status not in _STATUSES
                or not row.run_id.strip() or not row.dependency_group.strip()):
            raise ValueError(f"Invalid run or undeclared preparation: {row.run_id}")
        if row.status == "pass" and any(value is not None and not isfinite(value)
                                        for value in (row.measured_purity_pct, row.reference_purity_pct)):
            raise ValueError(f"Nonfinite purity for {row.run_id}")
        if row.run_id in by_run:
            if row != by_run[row.run_id]:
                raise ValueError(f"Conflicting rows for run {row.run_id}")
            continue
        by_run[row.run_id] = row
        repeats[row.preparation].append(row)
    for pair, covariance in covariances.items():
        groups = {by_run[run].dependency_group for run in pair if run in by_run}
        if len(groups) > 1 and covariance != 0:
            raise ValueError(f"Covariance contradicts independent groups for {sorted(pair)}")
    return repeats


def summarize_accuracy(
    rows: Iterable[AccuracyRow], *, preparations: Iterable[Preparation],
    exclusions: Mapping[Preparation, str] | None = None, tolerance_pp: float | None = None,
    k: float | None = None, error_covariances_pp2: Mapping[frozenset[str], float] | None = None,
) -> dict:
    """Score a frozen preparation roster, including entries with no runs.

    Summaries are keyed by (method, vendor, D1 stratum) and show the effective
    tolerance_pp. The top level shows the override and both decision 13 targets.
    Declare exclusions by preparation and give a reason. Identical copies
    of a run ID collapse. A passing run needs finite measured and reference
    purities to be scored.

    The uncertainty model assumes each measurement and its reference are
    independent for the variance of their difference. Reference u must already
    be standard uncertainty in percentage points, with its basis documented.
    Declare a dependency group on each run. Group IDs must capture all common
    error sources, including shared references. Count eligible groups from all
    runs of eligible preparations and result groups from runs used in the means.
    Supply total covariance of signed errors in pp² for each distinct run pair
    within a group, keyed by a frozenset of run IDs. This includes covariance
    from reused references and shared calibrants and may be explicitly zero.
    Distinct groups declare independence. The comparison uses |error|/(k*u)
    with propagated standard u and a caller-supplied k.
    """
    if tolerance_pp is not None and not _nonnegative(tolerance_pp):
        raise ValueError("Tolerance must be finite and nonnegative")
    exclusions, covariances = dict(exclusions or {}), error_covariances_pp2 or {}
    repeats = _roster_runs(rows, preparations, exclusions, k, covariances)

    results, strata, strata_repeats = [], defaultdict(list), defaultdict(list)
    for preparation, runs in repeats.items():
        usable, result = _preparation_result(preparation, runs, exclusions.get(preparation), k, covariances)
        key = preparation.method, preparation.vendor, preparation.d1_stratum
        results.append(result)
        strata[key].append(result)
        if usable:
            strata_repeats[key].append(usable)
    return {
        "estimator": ESTIMATOR, "tolerance_pp": tolerance_pp, "targets_pp": dict(_TARGETS_PP), "k": k,
        "uncertainty_comparison": "standard-u comparison with factor k" if k is not None else None,
        "preparations": results,
        "summaries": {
            key: _summarize(group, strata_repeats[key],
                            _TARGETS_PP[key[0]] if tolerance_pp is None else tolerance_pp, covariances)
            for key, group in strata.items()
        },
    }


def _independent_units(results):
    by_group = defaultdict(list)
    for index, result in enumerate(results):
        for group in result["dependency_groups"]:
            by_group[group].append(index)

    unseen = set(range(len(results)))
    units = []
    for index in range(len(results)):
        if index not in unseen:
            continue
        unseen.remove(index)
        pending, unit = [index], []
        while pending:
            result = results[pending.pop()]
            unit.append(result)
            for group in result["dependency_groups"]:
                for neighbor in by_group.pop(group, []):
                    if neighbor in unseen:
                        unseen.remove(neighbor)
                        pending.append(neighbor)
        units.append(unit)
    return units


def _gate(passed, reason, **numbers):
    status = "open" if passed is None else "pass" if passed else "fail"
    return {"status": status, "reason": reason, **numbers}


def _uncertainty_gate(scored, k, k_basis):
    reasons, missing, misses = [], [], []
    if k is None:
        reasons.append("Missing coverage factor k")
    if _missing(k_basis):
        reasons.append("Missing coverage basis k_basis")
    if not scored:
        reasons.append("No scored preparations")
    compared = 0
    for result in scored:
        standard_u = result["error_standard_u_pp"]
        preparation_id = result["preparation_id"]
        if standard_u is None:
            missing.append({"preparation_id": preparation_id, "reason": result["error_u_reason"]})
        elif k is not None:
            compared += 1
            if abs(result["error_pp"]) > k * standard_u:
                misses.append(preparation_id)
    if missing:
        reasons.append("Scored preparations lack standard uncertainty")
    agreeing = compared - len(misses)
    fraction = _fraction(agreeing, len(scored)) if compared == len(scored) else None
    return _gate(None if reasons else fraction >= 0.90,
                 ". ".join(reasons) if reasons else "At least 90% of scored preparations must agree within k*u",
                 k=k, k_basis=k_basis, scored_preparations=len(scored), compared_preparations=compared,
                 agreeing_preparations=agreeing, fraction=fraction, minimum_fraction=0.90,
                 misses=misses, missing_uncertainty=missing)


def _method_gates(results, target, k, k_basis):
    scored = [result for result in results if result["error_pp"] is not None]
    median_error = median(abs(result["error_pp"]) for result in scored) if scored else None
    units = _independent_units(results)
    scored_units = sum(any(result["error_pp"] is not None for result in unit) for unit in units)
    successes = sum(all(result["error_pp"] is not None and abs(result["error_pp"]) < target
                        and result["eligible_repeats"] == sum(result["statuses"].values())
                        for result in unit) for unit in units)
    n = len(units)
    p = sum(comb(n, count) for count in range(successes, n + 1)) / (1 << n) if n else None
    assessed = _fraction(len(scored), len(results))
    gates = {
        "accuracy": _gate(None if median_error is None else median_error <= target,
                          "No scored preparations" if median_error is None else
                          "Median absolute error must be at most the method target",
                          median_absolute_error_pp=median_error, target_pp=target, scored_preparations=len(scored)),
        "independent_units": _gate(True if scored_units >= 8 else None,
                                   "At least 8 independent units must have a scored preparation",
                                   scored_units=scored_units, eligible_units=n, minimum_scored_units=8),
        "sign_test": _gate(None if p is None else p <= 0.05,
                           "No eligible units" if p is None else "One-sided exact p must be at most 0.05",
                           successes=successes, failures=n - successes, units=n, p_value=p, alpha=0.05),
        "fraction_assessed": _gate(None if assessed is None else assessed >= 0.90,
                                  "No eligible preparations" if assessed is None else
                                  "At least 90% of eligible preparations must be scored",
                                  scored_preparations=len(scored), eligible_preparations=len(results),
                                  fraction=assessed, minimum_fraction=0.90),
        "uncertainty_agreement": _uncertainty_gate(scored, k, k_basis),
    }
    statuses = {gate["status"] for gate in gates.values()}
    status = "fail" if "fail" in statuses else "open" if "open" in statuses else "pass"
    return {"status": status, "gates": gates}


def evaluate_release_gates(
    rows: Iterable[AccuracyRow], *, preparations: Iterable[Preparation],
    exclusions: Mapping[Preparation, str] | None = None, k: float | None = None,
    k_basis: str | None = None, error_covariances_pp2: Mapping[frozenset[str], float] | None = None,
) -> dict:
    """Evaluate decision 13 on the frozen roster, returning a report keyed by method.

    Use the roster, dependency and covariance declarations of summarize_accuracy.
    Pool runs across vendors and D1 strata by (method, preparation_id).
    An exclusion on any stratum excludes the entire physical preparation.
    All remaining units enter the sign test, which requires every run in a
    successful unit to pass with finite measured and reference purities.
    An empty eligible roster leaves its gates open.
    Supply both a coverage factor k and its documented basis for uncertainty
    agreement. The report lists observed misses even when other evidence is missing.
    """
    if k_basis is not None and not isinstance(k_basis, str):
        raise ValueError("k_basis must be a string or None")
    exclusions, covariances = dict(exclusions or {}), error_covariances_pp2 or {}
    repeats = _roster_runs(rows, preparations, exclusions, k, covariances)
    excluded = {(prep.method, prep.preparation_id) for prep in exclusions}
    pooled = {}
    for preparation, runs in repeats.items():
        key = preparation.method, preparation.preparation_id
        if key not in excluded:
            pooled.setdefault(key, (preparation, []))[1].extend(runs)
    by_method = defaultdict(list)
    for preparation, runs in pooled.values():
        _, result = _preparation_result(preparation, runs, None, k, covariances)
        by_method[preparation.method].append(result)
    return {method: _method_gates(by_method[method], target, k, k_basis)
            for method, target in _TARGETS_PP.items()}
