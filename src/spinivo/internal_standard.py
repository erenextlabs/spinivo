"""Calculate same-tube qNMR purity and its GUM uncertainty budget."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field, replace
from fractions import Fraction
from itertools import combinations
from math import fsum, isfinite, sqrt
from numbers import Rational, Real
from typing import NamedTuple, TYPE_CHECKING

if TYPE_CHECKING:
    from .acquisition import CheckResult
    from .integral_components import RegionResult


@dataclass(frozen=True)
class Quantity:
    value: float | None
    standard_u: float | None
    unit: str
    source: str
    source_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class Correlation:
    inputs: tuple[str, str]
    r: float | None
    source: str
    independence_reason: str | None = None


@dataclass(frozen=True)
class Covariance:
    """Use the product of the input units, listed in input order."""

    inputs: tuple[str, str]
    value: float | None
    units: tuple[str, str]
    source: str
    independence_reason: str | None = None


@dataclass(frozen=True)
class CoverageFactor:
    value: float | None
    source: str


@dataclass(frozen=True)
class Prerequisite:
    name: str
    status: str
    reason: str

    @property
    def message(self) -> str:
        return f"{self.name}: {self.reason}"


@dataclass(frozen=True)
class InternalStandardAssay:
    """Use common units within each mass pair and mass fractions for purity.

    Area sources follow the order of the corresponding area values.
    Supply one aggregate correlation or covariance per pair sharing a source ID.
    Declare the signal names for each region when linking measured areas.
    """

    assay_id: str | None
    analyte_tube_id: str | None
    standard_tube_id: str | None
    analyte_areas: tuple[float | None, ...] | None
    standard_areas: tuple[float | None, ...] | None
    N_a: int | None
    N_s: int | None
    M_a: Quantity | None
    M_s: Quantity | None
    m_a: Quantity | None
    m_s: Quantity | None
    w_s: Quantity | None
    R_relative_u: float | None
    R_u_source: str
    correlations: tuple[Correlation | Covariance, ...] = ()
    k: CoverageFactor | None = None
    prerequisites: tuple[Prerequisite, ...] = ()
    analyte_area_sources: tuple[str | None, ...] | None = None
    standard_area_sources: tuple[str | None, ...] | None = None
    area_unit: str | None = None
    R_source_ids: tuple[str, ...] = ()
    analyte_signals: tuple[str, ...] = ()
    standard_signals: tuple[str, ...] = ()
    _linked_prerequisites: tuple[Prerequisite, ...] = field(default=(), repr=False, compare=False)


@dataclass(frozen=True)
class InputContribution:
    name: str
    quantity: Quantity
    sensitivity: float
    contribution: float
    share: float | None


@dataclass(frozen=True)
class CovarianceContribution:
    correlation: Correlation | Covariance
    term: float
    share: float | None
    r: float


@dataclass(frozen=True)
class UncertaintyBudget:
    """Shares are fractions of the total variance and undefined when it is zero."""

    inputs: tuple[InputContribution, ...]
    covariances: tuple[CovarianceContribution, ...]
    variance: float
    standard_u: float
    relative_standard_u: float
    k: CoverageFactor | None
    expanded_u: float | None
    expanded_u_reason: str | None


@dataclass(frozen=True)
class InternalStandardResult:
    assay_id: str | None
    status: str
    I_a: float | None
    I_s: float | None
    R: float | None
    value: float | None = None
    budget: UncertaintyBudget | None = None
    reasons: tuple[str, ...] = ()


class _InternalStandardCalculation(NamedTuple):
    I_a: float | None
    I_s: float | None
    R: float | None
    value: float | None
    budget: UncertaintyBudget | None
    reasons: tuple[str, ...]


CALCULATION_ONLY = "calculation only: no spectral evidence"
_QUANTITIES = ("M_a", "M_s", "m_s", "m_a", "w_s")
_FRACTION_UNITS = ("fraction", "1", "g/g", "mg/mg", "kg/kg")
_REQUIRED_CHECKS = ("recovery", "snr", "resolution")
BLANK_SIGNAL_REASON = "Region signal name is blank"


def link_spectrum_checks(assay: InternalStandardAssay, *, run: str,
                         analyte_regions: Iterable[tuple[str, RegionResult]],
                         standard_regions: Iterable[tuple[str, RegionResult]],
                         acquisition: Iterable[CheckResult]) -> InternalStandardAssay:
    """Attach nominal areas and the evidence that qualifies them for purity."""
    analyte_regions, standard_regions = tuple(analyte_regions), tuple(standard_regions)
    prerequisites = []
    signals = ()
    area_sources = {}
    region_units = {}
    for role, regions in (("analyte", analyte_regions), ("standard", standard_regions)):
        area_sources[role], expected = _link_regions(role, regions, getattr(assay, f"{role}_signals"),
                                                    run, prerequisites, region_units)
        signals += expected
    area_unit = _link_units(region_units, prerequisites)
    _link_checks("", signals, run, acquisition, prerequisites)
    return _replace_link_prerequisites(
        assay, prerequisites, analyte_areas=tuple(region.nominal_area for _, region in analyte_regions),
        standard_areas=tuple(region.nominal_area for _, region in standard_regions),
        analyte_area_sources=area_sources["analyte"], standard_area_sources=area_sources["standard"], area_unit=area_unit,
    )


def _replace_link_prerequisites(assay, prerequisites, **changes):
    # Keep caller-supplied entries, even when they equal a generated entry.
    previous = {id(check) for check in assay._linked_prerequisites}
    supplied = tuple(check for check in assay.prerequisites if id(check) not in previous)
    linked = tuple(prerequisites)
    return replace(assay, prerequisites=(*supplied, *linked), _linked_prerequisites=linked, **changes)


def _link_regions(role, regions, expected, run, prerequisites, region_units):
    if not isinstance(expected, (tuple, list)) or not expected or any(not _has_text(signal) for signal in expected):
        prerequisites.append(Prerequisite(f"{role} signals", "cannot be assessed", "Declare nonblank signal names"))
        expected = ()
    missing = list(expected)
    sources = []
    for index, (signal, region) in enumerate(regions, 1):
        name = f"{signal} {role} region {index}" if _has_text(signal) else f"{role} region {index}"
        if signal not in missing:
            reason = BLANK_SIGNAL_REASON if not _has_text(signal) else "Unexpected or repeated region"
            prerequisites.append(Prerequisite(name, "cannot be assessed", reason))
            sources.append(None)
            continue
        missing.remove(signal)
        region_run = _link_region(name, region, run, prerequisites, region_units)
        sources.append(f"{region_run} {name}" if region_run is not None else None)
    prerequisites.extend(Prerequisite(f"{signal} {role} region", "cannot be assessed", "Missing region")
                         for signal in missing)
    return tuple(sources), tuple(expected)


def _link_region(name, region, run, prerequisites, region_units):
    """Keep the region's own evidence and return the run its area is attributed to."""
    region_run = region.run
    if not _has_text(region_run):
        region_run = None
    elif region_run != run:
        prerequisites.append(Prerequisite(name, "cannot be assessed",
                                          f"Run mismatch: expected {run!r}, got {region_run!r}"))
        name = f"{region_run} {name}"
    reason = ". ".join(f"{check}: {cause}" for check, cause in (*region.flags, *region.gaps))
    prerequisites.append(Prerequisite(name, region.status, reason or "All region checks pass"))
    if region_run is None:
        prerequisites.append(Prerequisite(name, "cannot be assessed", "region run unknown"))
    if region_run is None or region_run == run:
        unit = region.area.unit
        region_units[name] = unit if _has_text(unit) else None
    return region_run


def _link_units(region_units, prerequisites):
    """Return the one nonblank area unit shared by every linked region."""
    area_units = set(region_units.values())
    if len(area_units) > 1 or None in area_units:
        units = ", ".join(f"{name} = {unit!r}" for name, unit in region_units.items())
        problem = "Region units are missing" if None in area_units else "Region units differ"
        prerequisites.append(Prerequisite("linked area units", "cannot be assessed", f"{problem} ({units})"))
        return None
    return area_units.pop() if area_units else None


def _link_checks(prefix, signals, run, acquisition, prerequisites):
    """Attach each contributing signal's checks from the linked run, and only those."""
    linked = set()
    for check in acquisition:
        if check.run == run and not _has_text(check.signal):
            prerequisites.append(Prerequisite(f"{prefix}{check.check}", "cannot be assessed", "Check signal name is blank"))
            continue
        if check.signal not in signals:
            continue
        name = f"{prefix}{check.signal} {check.check}"
        if check.run != run:
            prerequisites.append(Prerequisite(name, "cannot be assessed",
                                              f"Run mismatch: expected {run!r}, got {check.run!r}"))
            continue
        status, reason = check.status, check.reason
        if status == "pass" and check.check == "recovery":
            inputs = {item.name: item.source for item in check.inputs}
            angle_source = inputs.get("calibrated_flip_angle_deg")
            if _has_text(angle_source) and angle_source.endswith(" / calibrated_90_us"):
                angle_source = None  # A derived angle's formula is not its calibration source.
            calibration = inputs.get("calibrated_90_us", angle_source)
            missing = [f"{field} is missing" for field, source in
                       (("t1_source", inputs.get("t1_s")), ("calibration_source", calibration)) if not _has_text(source)]
            if missing:
                status, reason = "cannot be assessed", ". ".join(missing)
        prerequisites.append(Prerequisite(name, status, reason))
        linked.add((check.signal, check.check))
    prerequisites.extend(Prerequisite(f"{prefix}{signal} {check}", "cannot be assessed", f"Missing {check} check")
                         for signal in dict.fromkeys(signals) for check in _REQUIRED_CHECKS if (signal, check) not in linked)


def _has_text(value):
    return isinstance(value, str) and bool(value.strip())


def _finite(value):
    try:
        return isinstance(value, Real) and not isinstance(value, bool) and isfinite(value)
    except OverflowError:
        return False


def _validate_number(name, value, reasons, *, zero_ok=False):
    if _finite(value) and (value >= 0 if zero_ok else value > 0):
        return True
    bound = "nonnegative" if zero_ok else "positive"
    reasons.append(f"{name} is missing" if value is None else f"{name} must be finite and {bound}")
    return False


def _area_sum(name, areas, reasons):
    if not areas:
        reasons.append(f"{name} has no region areas")
        return None
    before = len(reasons)
    for index, area in enumerate(areas, 1):
        if not _finite(area):
            reasons.append(f"{name} region {index} area {'is missing' if area is None else 'must be finite'}")
    if len(reasons) != before:
        return None
    try:
        total = fsum(areas)
    except OverflowError:
        reasons.append(f"{name} area sum is outside the numeric range")
        return None
    if total <= 0:
        reasons.append(f"{name} area sum must be positive")
    return total


def _is_psd(names, correlations):
    zero, one = Fraction(0), Fraction(1)
    matrix = [[one if left == right else correlations.get(frozenset((left, right)), zero)
               for right in names] for left in names]
    while matrix:
        # Pivot on the largest diagonal before taking the Schur complement.
        pivot = max(range(len(matrix)), key=lambda i: matrix[i][i])
        diagonal = matrix[pivot][pivot]
        if diagonal <= 0:
            return not any(value for row in matrix for value in row)
        rest = [i for i in range(len(matrix)) if i != pivot]
        matrix = [[matrix[i][j] - matrix[i][pivot] * matrix[pivot][j] / diagonal for j in rest] for i in rest]
    return True


def _fraction(value):
    return (Fraction(int(value.numerator), int(value.denominator)) if isinstance(value, Rational)
            else Fraction(str(float(value))))


def _ratio_quantity(R, assay):
    u = abs(R) * assay.R_relative_u if R is not None and _finite(assay.R_relative_u) else None
    return Quantity(R, u, "ratio", assay.R_u_source, assay.R_source_ids)


def _validate_correlations(correlations, quantities, reasons):
    if not isinstance(correlations, (tuple, list)) or any(not isinstance(c, (Correlation, Covariance)) for c in correlations):
        reasons.append("correlations must contain Correlation records or Covariance records")
        return ()
    before = len(reasons)
    entries, seen, groups = {}, set(), {}
    for name, quantity in quantities.items():
        if quantity is None:
            continue
        ids = quantity.source_ids
        if not isinstance(ids, (tuple, list)) or any(not _has_text(i) or i.strip() != i for i in ids):
            reasons.append(f"{name}.source_ids must contain nonblank source identifiers without outer whitespace")
            continue
        for source_id in dict.fromkeys(ids):
            groups.setdefault(source_id, []).append(name)
    shared = {frozenset(pair): source_id for source_id, names in groups.items() for pair in combinations(names, 2)}
    for record in correlations:
        inputs = record.inputs
        raw = isinstance(record, Covariance)
        kind = "covariance" if raw else "correlation"
        if not isinstance(inputs, (tuple, list)) or any(not isinstance(item, str) for item in inputs):
            reasons.append(f"{kind}.inputs must contain two input names as strings")
            continue
        pair = frozenset(inputs)
        name = f"{kind} {', '.join(inputs)}"
        if len(inputs) != 2 or len(pair) != 2:
            reasons.append(f"{name} must name two distinct inputs")
        elif pair.difference(quantities):
            reasons.append(f"{name} contains an unknown budget input")
        elif pair in seen:
            reasons.append(f"{name} is a duplicate pair")
        else:
            seen.add(pair)
            r = _covariance_r(record, quantities, name, reasons) if raw else record.r
            if not raw and (not _finite(r) or abs(r) > 1):
                reasons.append(f"{name} requires a finite r in [-1, 1]")
                r = None
            if r is not None:
                if r == 0 and pair in shared and not _has_text(record.independence_reason):
                    reasons.append(f"{name} requires a justification of independence")
                elif r != 0 and record.independence_reason is not None:
                    reasons.append(f"{name} cannot claim independence with a nonzero value")
                entries[pair] = (record, _fraction(r))
        if not _has_text(record.source):
            reasons.append(f"{name} source is missing")
    names = tuple(quantities)
    for pair, source_id in shared.items():
        if pair not in seen:
            joined = ", ".join(name for name in names if name in pair)
            reasons.append(f"missing covariance for {joined} (shared source {source_id})")
    coefficients = {pair: r for pair, (_, r) in entries.items()}
    if len(reasons) == before and len(entries) == len(seen) and not _is_psd(names, coefficients):
        conflict = names
        for name in names:
            subset = tuple(item for item in conflict if item != name)
            if not _is_psd(subset, coefficients):
                conflict = subset
        reasons.append(f"correlation matrix is not positive semidefinite for {', '.join(conflict)}")
    return tuple(entries.values())


def _covariance_r(record, quantities, name, reasons):
    left, right = (quantities[key] for key in record.inputs)
    if left is None or right is None:
        return None
    if not isinstance(record.units, (tuple, list)) or tuple(record.units) != (left.unit, right.unit):
        reasons.append(f"{name} units must match the input units in order")
        return None
    if not _finite(record.value):
        reasons.append(f"{name} requires a finite value")
        return None
    if not all(_finite(q.standard_u) and q.standard_u >= 0 for q in (left, right)):
        return None
    limit = _fraction(left.standard_u) * _fraction(right.standard_u)
    value = _fraction(record.value)
    if abs(value) > limit:
        reasons.append(f"{name} exceeds the product of the standard uncertainties")
        return None
    return value / limit if limit else Fraction(0)


def _uncertainty_budget(value, quantities, correlations, k, negative_inputs=("M_s", "m_a")):
    sensitivities = {name: (-value if name in negative_inputs else value) / quantity.value
                     for name, quantity in quantities.items()}
    contributions = {name: sensitivities[name] * quantity.standard_u for name, quantity in quantities.items()}
    squares = {name: contribution * contribution for name, contribution in contributions.items()}
    terms = []
    for correlation, r in correlations:
        left, right = correlation.inputs
        terms.append(2 * (contributions[left] * contributions[right]) * float(r))
    variance = fsum((*squares.values(), *terms))
    if not isfinite(variance) or variance < 0:
        raise ValueError("combined variance must be finite and nonnegative")
    standard_u = sqrt(variance)
    relative_u = standard_u / value
    expanded_u = k.value * standard_u if k is not None else None
    if not isfinite(relative_u) or (expanded_u is not None and not isfinite(expanded_u)):
        raise ValueError("relative or expanded uncertainty is outside the numeric range")
    inputs = tuple(InputContribution(name, quantity, sensitivities[name], contributions[name],
                                     squares[name] / variance if variance else None)
                   for name, quantity in quantities.items())
    covariances = tuple(CovarianceContribution(record, term, term / variance if variance else None, float(r))
                        for (record, r), term in zip(correlations, terms))
    return UncertaintyBudget(inputs, covariances, variance, standard_u, relative_u, k, expanded_u,
                             "coverage factor k was not supplied" if k is None else None)


def calculate_internal_standard(assay: InternalStandardAssay) -> InternalStandardResult:
    """Calculate an entered or edited record; only spinivo.pipeline.assess_internal_standard can pass."""
    I_a, I_s, R, value, budget, reasons = _calculate_internal_standard(assay)
    return InternalStandardResult(assay.assay_id, "cannot be assessed", I_a, I_s, R, value, budget,
                                  reasons or (CALCULATION_ONLY,))


def _calculate_internal_standard(assay: InternalStandardAssay) -> _InternalStandardCalculation:
    """Preserve computable areas and their ratio even when purity is blocked."""
    reasons = []
    I_a = _area_sum("I_a (analyte_areas)", assay.analyte_areas, reasons)
    I_s = _area_sum("I_s (standard_areas)", assay.standard_areas, reasons)
    R = I_a / I_s if I_a is not None and I_s else None
    if R is not None and not isfinite(R):
        reasons.append("R is outside the numeric range")
        R = None
    elif R == 0 and I_a > 0 and I_s > 0:
        reasons.append("R is too small to represent")

    for role in ("analyte", "standard"):
        name = f"{role}_area_sources"
        sources = getattr(assay, name)
        if sources is None:
            continue
        if not isinstance(sources, (tuple, list)):
            reasons.append(f"{name} must be a tuple or list of area sources")
            continue
        area_count = len(getattr(assay, f"{role}_areas") or ())
        source_count = len(sources)
        if source_count != area_count:
            entries = "entry" if source_count == 1 else "entries"
            areas = "area" if area_count == 1 else "areas"
            reasons.append(f"{name} has {source_count} {entries} for {area_count} {areas}")

    for name in ("assay_id", "analyte_tube_id", "standard_tube_id", "R_u_source"):
        if not _has_text(getattr(assay, name)):
            reasons.append(f"{name} is missing")
    tube_a, tube_s = assay.analyte_tube_id, assay.standard_tube_id
    if _has_text(tube_a) and _has_text(tube_s) and tube_a != tube_s:
        reasons.append("analyte_tube_id and standard_tube_id differ")
    for name in ("N_a", "N_s"):
        count = getattr(assay, name)
        if _validate_number(name, count, reasons) and count % 1:
            reasons.append(f"{name} must be a whole count")
    quantities = {name: getattr(assay, name) for name in _QUANTITIES}
    for name, quantity in quantities.items():
        if quantity is None:
            reasons.append(f"{name} is missing")
            continue
        if _validate_number(f"{name}.value", quantity.value, reasons) and name == "w_s" and quantity.value > 1:
            reasons.append("w_s.value must be at most 1")
        _validate_number(f"{name}.standard_u", quantity.standard_u, reasons, zero_ok=True)
        if name == "w_s":
            if not isinstance(quantity.unit, str) or quantity.unit not in _FRACTION_UNITS:
                reasons.append(f"w_s.unit must be one of {', '.join(_FRACTION_UNITS)}")
        elif not _has_text(quantity.unit):
            reasons.append(f"{name}.unit is missing")
        if not _has_text(quantity.source):
            reasons.append(f"{name}.source is missing")
    for left, right in (("M_a", "M_s"), ("m_a", "m_s")):
        a, b = quantities[left], quantities[right]
        if a and b and _has_text(a.unit) and _has_text(b.unit) and a.unit != b.unit:
            reasons.append(f"{left} and {right} must use the same unit")
    _validate_number("R_relative_u", assay.R_relative_u, reasons, zero_ok=True)
    quantities = {"R": _ratio_quantity(R, assay), **quantities}
    correlations = _validate_correlations(assay.correlations, quantities, reasons)
    if assay.k is not None:
        _validate_number("k.value", assay.k.value, reasons)
        if not _has_text(assay.k.source):
            reasons.append("k.source is missing")
    reasons.extend(check.message for check in assay.prerequisites if check.status != "pass")

    if not reasons:
        try:
            value = (R * (assay.N_s / assay.N_a) * (assay.M_a.value / assay.M_s.value)
                     * (assay.m_s.value / assay.m_a.value) * assay.w_s.value)
            in_range = isfinite(value) and value > 0
        except ArithmeticError:
            in_range = False
        if in_range:
            try:
                budget = _uncertainty_budget(value, quantities, correlations, assay.k)
            except (ArithmeticError, ValueError) as error:
                reasons.append(f"uncertainty budget: {error}")
            else:
                return _InternalStandardCalculation(I_a, I_s, R, value, budget, ())
        else:
            reasons.append("calculated purity is outside the numeric range (R, N_s, N_a, M_a, M_s, m_s, m_a, w_s)")
    return _InternalStandardCalculation(I_a, I_s, R, None, None, tuple(reasons))
