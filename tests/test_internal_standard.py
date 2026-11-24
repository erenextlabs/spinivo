import json
import unittest
from dataclasses import replace
from fractions import Fraction
from functools import lru_cache
from math import inf, isclose, nan, nextafter
from pathlib import Path

from spinivo.acquisition import RecoveryEvidence, SignalEvidence, SignalMeasurement, check_acquisition
from spinivo.integral_components import ComponentEvidence, Line, check_region
from spinivo.internal_standard import (
    CALCULATION_ONLY, Correlation, Covariance, CoverageFactor, InternalStandardAssay, Prerequisite, Quantity,
    _calculate_internal_standard, calculate_internal_standard, link_spectrum_checks,
)
from spinivo.processing import process
from tests.integral_fixtures import NOISE_BOUNDS_HZ, SETTINGS, load_cases, make_raw, ppm_bounds


QUANTITY_NAMES = ("M_a", "M_s", "m_a", "m_s", "w_s")
NOISE = ppm_bounds(NOISE_BOUNDS_HZ)
NOISE_GAP = "snr: Missing noise window. neighbour: Missing noise window"


@lru_cache(maxsize=None)
def clean_regions():
    """Return the clean fixture raw and its region checked with and without the noise window."""
    case = next(case for case in load_cases()["cases"] if case["id"] == "clean_singlet")
    raw, bounds = make_raw(case), ppm_bounds(case["bounds"])
    spectrum = process(raw, SETTINGS, regions_ppm=(bounds,))
    evidence = ComponentEvidence(target_lines=(Line(5.0, 1.0, 2.0),), target_source="fixture",
                                 satellite_convention="absent", satellite_source="fixture")
    return raw, *(check_region(spectrum, bounds, evidence, raw=raw, noise_bounds_ppm=noise) for noise in (NOISE, None))


def acquisition_checks(raw, signals, *, noise=NOISE, no_t1_source=(), calibration_source="pulse log", narrow=(),
                       via_pulse_width=False):
    """Run the three acquisition checks on one run; narrow lines flag resolution."""
    raw = replace(raw, relaxation_delay_s=60.0, pulse_width_us=10.0)
    evidence = RecoveryEvidence(calibrated_flip_angle_deg=None if via_pulse_width else 90.0,
                                calibrated_90_us=10.0 if via_pulse_width else None, calibration_source=calibration_source)
    return check_acquisition(raw, [SignalEvidence(signal, 2.0, None if signal in no_t1_source else "T1 log",
                                                  SignalMeasurement(3000.0, 1.0 if noise else None,
                                                                    0.5 if signal in narrow else 2.0, "fixture", noise))
                                   for signal in signals], evidence)


def load_assay_case():
    path = Path(__file__).parent / "fixtures" / "internal_standard_cases.json"
    case = json.loads(path.read_text())
    inputs = case["assay"]
    for name in QUANTITY_NAMES:
        fields = inputs[name]
        inputs[name] = Quantity(**{**fields, "source_ids": tuple(fields.get("source_ids", ()))})
    for name in ("analyte_areas", "standard_areas", "analyte_signals", "standard_signals", "R_source_ids"):
        inputs[name] = tuple(inputs.get(name, ()))
    inputs["correlations"] = tuple(Correlation(tuple(c["inputs"]), c["r"], c["source"]) for c in inputs["correlations"])
    inputs["k"] = CoverageFactor(**inputs["k"])
    return InternalStandardAssay(**inputs), case["expected"]


class InternalStandardTests(unittest.TestCase):
    def setUp(self):
        self.assay, self.expected = load_assay_case()
        self.original = calculate_internal_standard(self.assay)

    def assertRelative(self, actual, expected):
        self.assertTrue(isclose(actual, expected, rel_tol=1e-9), (actual, expected))

    def assertBlocked(self, assay, reason):
        result = calculate_internal_standard(assay)
        self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
        self.assertTrue(any(reason in text for text in result.reasons), result.reasons)
        return result

    def test_calculation_core_returns_values_without_a_status(self):
        calculated = _calculate_internal_standard(self.assay)
        self.assertFalse(hasattr(calculated, "status"))
        self.assertEqual((calculated.value, calculated.budget, calculated.reasons),
                         (self.original.value, self.original.budget, ()))
        blocked = _calculate_internal_standard(replace(self.assay, m_a=None))
        self.assertEqual((blocked.value, blocked.budget, blocked.reasons), (None, None, ("m_a is missing",)))

    def test_hand_worked_budget(self):
        result, expected = self.original, self.expected
        self.assertEqual((result.status, result.assay_id, result.reasons),
                         ("cannot be assessed", self.assay.assay_id, (CALCULATION_ONLY,)))
        for name in ("I_a", "I_s", "R"):
            self.assertRelative(getattr(result, name), expected[name])
        self.assertEqual(round(result.value, 10), expected["purity"])
        budget = result.budget
        for name in ("variance", "standard_u", "relative_standard_u", "expanded_u"):
            self.assertRelative(getattr(budget, name), expected[name])
        self.assertEqual(budget.k, self.assay.k)
        self.assertIsNone(budget.expanded_u_reason)
        self.assertEqual({item.name for item in budget.inputs}, set(expected["sensitivities"]))
        for item in budget.inputs:
            with self.subTest(input=item.name):
                self.assertIsInstance(item.quantity.source_ids, tuple)
                self.assertRelative(item.sensitivity, expected["sensitivities"][item.name])
                contribution = expected["contributions"][item.name]
                self.assertRelative(item.contribution, contribution)
                self.assertRelative(item.share, contribution ** 2 / expected["variance"])
                quantity = (Quantity(result.R, result.R * self.assay.R_relative_u, "ratio", self.assay.R_u_source)
                            if item.name == "R" else getattr(self.assay, item.name))
                self.assertEqual(item.quantity, quantity)
        covariance, = budget.covariances
        self.assertEqual(covariance.correlation, self.assay.correlations[0])
        self.assertRelative(covariance.term, expected["covariance_term"])
        self.assertRelative(covariance.share, expected["covariance_share"])

    def test_missing_facts_and_restoration(self):
        for name in ("assay_id", "analyte_tube_id", "standard_tube_id", "analyte_areas", "standard_areas",
                     "N_a", "N_s", "R_relative_u", "R_u_source", *QUANTITY_NAMES):
            with self.subTest(input=name):
                missing = replace(self.assay, **{name: None})
                self.assertBlocked(missing, name)
                restored = replace(missing, **{name: getattr(self.assay, name)})
                self.assertEqual(calculate_internal_standard(restored), self.original)
        for name in QUANTITY_NAMES:
            quantity = getattr(self.assay, name)
            for field in ("value", "standard_u", "source", "unit"):
                with self.subTest(input=name, field=field):
                    missing = replace(self.assay, **{name: replace(quantity, **{field: None})})
                    self.assertBlocked(missing, f"{name}.{field}")
                    self.assertEqual(calculate_internal_standard(replace(missing, **{name: quantity})), self.original)
        for field in ("value", "source"):
            with self.subTest(input="k", field=field):
                missing = replace(self.assay, k=replace(self.assay.k, **{field: None}))
                self.assertBlocked(missing, f"k.{field}")
                self.assertEqual(calculate_internal_standard(replace(missing, k=self.assay.k)), self.original)
        correlation, = self.assay.correlations
        for field, reason in (("r", "requires a finite r"), ("source", "source is missing")):
            with self.subTest(input="correlation", field=field):
                missing = replace(self.assay, correlations=(replace(correlation, **{field: None}),))
                self.assertBlocked(missing, f"correlation m_s, m_a {reason}")
                restored = replace(missing, correlations=self.assay.correlations)
                self.assertEqual(calculate_internal_standard(restored), self.original)

    def test_blockers_preserve_areas_and_ratio(self):
        correlation, = self.assay.correlations
        for changes, reason in (
            ({"standard_tube_id": "another tube"}, "analyte_tube_id and standard_tube_id differ"),
            ({"correlations": (replace(correlation, r=1.1),)}, "correlation m_s, m_a"),
            ({"correlations": (replace(correlation, inputs=("R", "N_s")),)}, "correlation R, N_s contains an unknown"),
            ({"correlations": (correlation, replace(correlation, inputs=("m_a", "m_s")))}, "duplicate pair"),
            ({"prerequisites": (Prerequisite("baseline", "flag", "distorted integral"),)}, "baseline: distorted integral"),
            ({"w_s": replace(self.assay.w_s, value=1.001)}, "w_s.value"),
            ({"w_s": replace(self.assay.w_s, unit="%")}, "w_s.unit"),
            ({"w_s": replace(self.assay.w_s, unit="ppm")}, "w_s.unit"),
            ({"M_a": replace(self.assay.M_a, value=inf)}, "M_a.value"),
            ({"M_s": replace(self.assay.M_s, standard_u=-0.01)}, "M_s.standard_u"),
            ({"N_a": 0}, "N_a"),
            ({"N_a": 10**400}, "N_a"),
            ({"R_relative_u": nan}, "R_relative_u"),
            ({"analyte_area_sources": ("H5", "I5", "extra")}, "analyte_area_sources has 3 entries for 2 areas"),
            ({"standard_area_sources": ()}, "standard_area_sources has 0 entries for 1 area"),
        ):
            with self.subTest(reason=reason):
                result = self.assertBlocked(replace(self.assay, **changes), reason)
                self.assertEqual((result.I_a, result.I_s, result.R),
                                 (self.original.I_a, self.original.I_s, self.original.R))

    def test_purity_overflow_preserves_areas_and_ratio(self):
        for denominator in (Fraction(1, 10**400), 1e-310):
            with self.subTest(denominator=denominator):
                assay = replace(self.assay, M_a=replace(self.assay.M_a, value=1),
                                M_s=replace(self.assay.M_s, value=denominator))
                result = self.assertBlocked(assay, "calculated purity is outside the numeric range")
                self.assertTrue(any("M_a, M_s" in reason for reason in result.reasons))
                self.assertEqual((result.I_a, result.I_s, result.R),
                                 (self.original.I_a, self.original.I_s, self.original.R))

    def test_correlation_matrix(self):
        weighing, = self.assay.correlations
        pairs = (("M_a", "m_s"), ("m_s", "w_s"), ("M_a", "w_s"))
        valid = tuple(Correlation(pair, 1, "fixture") for pair in pairs)
        passing = replace(self.assay, correlations=(*valid, replace(weighing, r=0, independence_reason="Separate tares")))
        self.assertEqual(calculate_internal_standard(passing).reasons, (CALCULATION_ONLY,))
        unrelated = Correlation(("M_s", "m_a"), 0.3, "fixture")
        singular = (replace(valid[0], r=0.6), replace(valid[1], r=0.8))
        self.assertEqual(calculate_internal_standard(replace(passing, correlations=(*singular, passing.correlations[-1]))).reasons,
                         (CALCULATION_ONLY,))
        for invalid in ((*valid[:2], replace(valid[2], r=-1)), valid[:2],
                        (*valid[:2], replace(valid[2], r=nextafter(1.0, 0.0)))):
            with self.subTest(correlations=invalid):
                blocked = replace(self.assay, correlations=(*invalid, unrelated, passing.correlations[-1]))
                self.assertBlocked(blocked, "not positive semidefinite for M_a, m_s, w_s")

    def test_required_covariance_evidence(self):
        missing = replace(self.assay, correlations=())
        self.assertBlocked(missing, "missing covariance for m_s, m_a (shared source balance:fixture-FID06)")
        self.assertEqual(calculate_internal_standard(replace(missing, correlations=self.assay.correlations)), self.original)
        weighing, = self.assay.correlations
        raw = Covariance(("m_a", "m_s"), 4.8e-6, ("mg", "mg"), weighing.source)
        for entry, reason in ((replace(raw, units=("g", "mg")), "units must match"),
                              (replace(raw, value=1.7e-5), "exceeds the product"),
                              (replace(weighing, r=0), "requires a justification")):
            self.assertBlocked(replace(missing, correlations=(entry,)), reason)
        self.assertBlocked(replace(missing, correlations=(weighing, raw)), "duplicate pair")
        self.assertEqual(calculate_internal_standard(replace(missing, correlations=(raw,))).budget.variance,
                         self.original.budget.variance)
        zero = replace(weighing, r=0, independence_reason="Independent tare tests")
        result = calculate_internal_standard(replace(missing, correlations=(zero,)))
        self.assertEqual(result.reasons, (CALCULATION_ONLY,))
        self.assertRelative(result.budget.variance, self.expected["variance"] - self.expected["covariance_term"])

    def test_correlation_input_names(self):
        correlation, = self.assay.correlations
        for inputs in (None, "m_a", ("m_s", None), ("m_s", 1), ([], "m_a")):
            with self.subTest(inputs=inputs):
                self.assertBlocked(replace(self.assay, correlations=(replace(correlation, inputs=inputs),)),
                                   "correlation.inputs")

    def test_fraction_units(self):
        for unit in ("fraction", "1", "g/g", "mg/mg", "kg/kg"):
            with self.subTest(unit=unit):
                assay = replace(self.assay, w_s=replace(self.assay.w_s, unit=unit))
                self.assertEqual(calculate_internal_standard(assay).value, self.original.value)

    def test_area_failures(self):
        for areas, total in (((), None), ((1.0, None), None), ((inf,), None), ((-1.0, 1.0), 0.0), ((-1.0,), -1.0)):
            with self.subTest(areas=areas):
                result = self.assertBlocked(replace(self.assay, analyte_areas=areas), "I_a")
                self.assertEqual((result.I_a, result.I_s), (total, self.original.I_s))
                self.assertEqual(result.R, total / self.original.I_s if total is not None else None)
        result = self.assertBlocked(replace(self.assay, standard_areas=(0.0,)), "I_s")
        self.assertEqual((result.I_a, result.I_s, result.R), (self.original.I_a, 0.0, None))

    def test_optional_coverage_factor(self):
        result = calculate_internal_standard(replace(self.assay, k=None))
        self.assertEqual((result.reasons, result.value), ((CALCULATION_ONLY,), self.original.value))
        self.assertEqual(result.budget.standard_u, self.original.budget.standard_u)
        self.assertIsNone(result.budget.expanded_u)
        self.assertIn("coverage factor k", result.budget.expanded_u_reason)


class LinkedEvidenceTests(unittest.TestCase):
    run_id, narrow = "FID06", ()
    signals = ("H5", "I5", "G5")
    names = ("H5 analyte region 1", "I5 analyte region 2", "G5 standard region 1")

    @classmethod
    def setUpClass(cls):
        cls.assay, _ = load_assay_case()
        cls.original = calculate_internal_standard(cls.assay)
        raw, noisy, quiet = clean_regions()
        cls.raw = replace(raw, source=Path(cls.run_id), pulse_program="single_pulse")
        cls.noisy, cls.quiet = (replace(region, run=cls.run_id) for region in (noisy, quiet))

    def regions(self, region=None, **changes):
        region = replace(region or self.noisy, **changes)
        areas = (*self.assay.analyte_areas, *self.assay.standard_areas)
        pairs = tuple((signal, replace(region, nominal_area=area)) for signal, area in zip(self.signals, areas))
        return {"analyte_regions": pairs[:2], "standard_regions": pairs[2:]}

    def checks(self, **changes):
        return {"acquisition": acquisition_checks(self.raw, self.signals, narrow=self.narrow, **changes)}

    def link(self, starting=None, **changes):
        if starting is None:
            starting = replace(self.assay, analyte_areas=None, standard_areas=None)
        return link_spectrum_checks(starting, run=self.run_id, **{**self.regions(), **self.checks(), **changes})

    def test_region_without_noise_keeps_area_sources(self):
        raw, _, quiet = clean_regions()
        run = str(raw.source)
        linked = link_spectrum_checks(self.assay, run=run, **self.regions(quiet),
                                      acquisition=acquisition_checks(raw, self.signals, noise=None))
        self.assertIsNone(quiet.snr.details)
        self.assertEqual((*linked.analyte_area_sources, *linked.standard_area_sources),
                         tuple(f"{run} {name}" for name in self.names))
        self.assertFalse(any(check.reason == "region run unknown" for check in linked.prerequisites))
        self.assertIn(f"{self.names[0]}: {NOISE_GAP}", calculate_internal_standard(linked).reasons)

    def test_required_regions_and_checks(self):
        for key, regions in self.regions().items():
            for index, (signal, _) in enumerate(regions):
                role = key.removesuffix("_regions")
                with self.subTest(role=role, signal=signal):
                    linked = self.link(**{key: regions[:index] + regions[index + 1:]})
                    self.assertLinkBlocked(linked, f"{signal} {role} region: Missing region")
        checks = self.checks()["acquisition"]
        for index, check in enumerate(checks):
            with self.subTest(signal=check.signal, check=check.check):
                linked = self.link(acquisition=checks[:index] + checks[index + 1:])
                self.assertLinkBlocked(linked, f"{check.signal} {check.check}: Missing {check.check} check")

    def assertLinkBlocked(self, assay, reason):
        result = calculate_internal_standard(assay)
        self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
        self.assertIn(reason, result.reasons)

    def test_blank_signal_names_block(self):
        analyte = self.regions()["analyte_regions"]
        checks = tuple(replace(check, signal=" ") if check.signal == "H5" else check
                       for check in self.checks()["acquisition"])
        linked = self.link(analyte_regions=((" ", analyte[0][1]), analyte[1]), acquisition=checks)
        self.assertLinkBlocked(linked, "analyte region 1: Region signal name is blank")
        self.assertLinkBlocked(linked, "recovery: Check signal name is blank")
        self.assertLinkBlocked(self.link(replace(self.assay, analyte_signals=(" ",))),
                               "analyte signals: Declare nonblank signal names")

    def test_foreign_flag_stays_with_its_run(self):
        foreign = replace(self.noisy, status="flag", flags=(("baseline", "distorted"),), run="other run")
        linked = self.link(**self.regions(foreign))
        self.assertLinkBlocked(linked, "H5 analyte region 1: Run mismatch: expected 'FID06', got 'other run'")
        self.assertEqual([p for p in linked.prerequisites if p.name == self.names[0]],
                         [Prerequisite(self.names[0], "cannot be assessed", "Run mismatch: expected 'FID06', got 'other run'")])
        self.assertIn(Prerequisite(f"other run {self.names[0]}", "flag", "baseline: distorted"), linked.prerequisites)
        self.assertEqual(linked.analyte_area_sources[0], f"other run {self.names[0]}")

    def test_relink_replaces_generated_prerequisites(self):
        supplied = Prerequisite("same tube", "pass", "Documented")
        starting = replace(self.assay, prerequisites=(supplied,))
        checks = self.checks()["acquisition"]
        blocked = self.link(starting, acquisition=(replace(checks[0], inputs=()), *checks[1:]))
        self.assertLinkBlocked(blocked, "H5 recovery: t1_source is missing. calibration_source is missing")
        corrected = self.link(blocked)
        self.assertEqual(corrected, self.link(starting))
        self.assertEqual(corrected.prerequisites[0], supplied)
        self.assertEqual(calculate_internal_standard(corrected), self.original)

    def test_recovery_pass_requires_each_source(self):
        for via_width, calibration in ((False, "calibrated_flip_angle_deg"), (True, "calibrated_90_us")):
            recovery, *other = self.checks(via_pulse_width=via_width)["acquisition"]
            self.assertEqual(calculate_internal_standard(self.link(acquisition=(recovery, *other))), self.original)
            for name, reason in (("t1_s", "t1_source"), (calibration, "calibration_source")):
                with self.subTest(calibration=calibration, missing=name):
                    inputs = tuple(item for item in recovery.inputs if item.name != name)
                    broken = self.link(acquisition=(replace(recovery, inputs=inputs), *other))
                    self.assertLinkBlocked(broken, f"H5 recovery: {reason} is missing")

    def test_missing_facts_block_only_dependent_results(self):
        units = "Region units are missing (" + ", ".join(f"{name} = None" for name in self.names) + ")"
        for self.narrow in ((), ("G5",)):
            baseline = self.link()
            expected = calculate_internal_standard(baseline)
            if self.narrow:
                self.assertIn("G5 resolution: Outside the pilot limit", expected.reasons)
            else:
                self.assertEqual(expected, self.original)
            cases = (
                ({"acquisition": self.checks()["acquisition"][1:]}, [("H5 recovery", "Missing recovery check")]),
                ({**self.regions(self.quiet), **self.checks(noise=None)},
                 [(name, NOISE_GAP) for name in self.names]
                 + [(f"{s} snr", "noise_sd is missing") for s in self.signals]),
                (self.checks(no_t1_source=("H5",)), [("H5 recovery", "t1_source is missing")]),
                (self.checks(calibration_source=" "), [(f"{s} recovery", "calibration_source is missing")
                                                       for s in self.signals]),
                *((self.regions(area=replace(self.noisy.area, unit=unit)), [("linked area units", units)])
                  for unit in ("", None)),
                (self.regions(run=None), [(name, "region run unknown") for name in self.names]),
            )
            for changes, added in cases:
                with self.subTest(narrow=self.narrow, fact=added[0]):
                    linked = self.link(**changes)
                    result = calculate_internal_standard(linked)
                    self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
                    added = [Prerequisite(name, "cannot be assessed", reason) for name, reason in added]
                    self.assertEqual([p for p in linked.prerequisites if p not in baseline.prerequisites], added)
                    self.assertIn(f"{added[0].name}: {added[0].reason}", result.reasons)
                    names = {p.name for p in added}
                    self.assertLessEqual({p for p in baseline.prerequisites if p.name not in names},
                                         set(linked.prerequisites))
                    for name in ("analyte_areas", "standard_areas"):
                        self.assertEqual(getattr(linked, name), getattr(baseline, name))
                    for role in ("analyte", "standard"):
                        name = f"{role}_area_sources"
                        sources = getattr(baseline, name)
                        if any(p.reason == "region run unknown" for p in added):
                            sources = (None,) * len(sources)
                        self.assertEqual(getattr(linked, name), sources)
                    self.assertEqual((result.I_a, result.I_s, result.R),
                                     (self.original.I_a, self.original.I_s, self.original.R))
                    self.assertEqual(calculate_internal_standard(self.link()), expected)
        self.assertEqual(baseline.analyte_area_sources, ("FID06 H5 analyte region 1", "FID06 I5 analyte region 2"))
