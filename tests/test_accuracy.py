import unittest
from dataclasses import replace
from math import sqrt

from spinivo import evaluate_release_gates
from spinivo.accuracy import ESTIMATOR, AccuracyRow, Preparation, summarize_accuracy


STRATUM = ("external_calibration", "JEOL", "long")
INTERNAL_STRATUM = ("internal_standard", "JEOL", "long")


def preparation(name, stratum=STRATUM):
    return Preparation(name, *stratum)


def result(prep, measured, name=None, **changes):
    row = AccuracyRow(name or prep.preparation_id, prep, "pass", prep.preparation_id, measured, 99,
                      measurement_u_pp=1, reference_u_pp=0.5, reference_coverage_basis="standard")
    return replace(row, **changes)


class AccuracyTests(unittest.TestCase):
    def test_preparation_weighting_signed_units_and_repeat_sd(self):
        a, b, c = (preparation(name) for name in "abc")
        rows = [result(a, 99, "a1"), result(a, 101, "a2"), result(b, 98.5), result(c, 97)]
        report = summarize_accuracy(rows, preparations=[a, b, c])
        summary = report["summaries"][STRATUM]
        self.assertEqual(report["estimator"], ESTIMATOR)
        self.assertIsNone(report["tolerance_pp"])
        self.assertEqual(report["preparations"][0]["measured_purity_pct"], 100)
        self.assertEqual(report["preparations"][0]["error_pp"], 1)
        self.assertAlmostEqual(report["preparations"][0]["repeat_sd_pp"], sqrt(2))
        self.assertIsNone(report["preparations"][1]["repeat_sd_pp"])
        self.assertEqual(summary["bias_pp"], -0.5)
        self.assertEqual(summary["median_absolute_error_pp"], 1)
        self.assertEqual(summary["agreement_fraction"], 1 / 3)
        strict = summarize_accuracy(rows, preparations=[a, b, c], tolerance_pp=0.25)
        self.assertEqual(strict["summaries"][STRATUM]["agreement_fraction"], 0)

    def test_method_targets_and_common_override(self):
        preps = [preparation("external"), preparation("internal", INTERNAL_STRATUM)]
        rows = [result(prep, 99.625) for prep in preps]
        report = summarize_accuracy(rows, preparations=preps)
        default = report["summaries"]
        self.assertEqual(report["targets_pp"], {"internal_standard": 0.50, "external_calibration": 0.75})
        self.assertEqual((default[STRATUM]["tolerance_pp"], default[INTERNAL_STRATUM]["tolerance_pp"]), (0.75, 0.50))
        self.assertEqual(default[STRATUM]["agreement_fraction"], 1)
        self.assertEqual(default[INTERNAL_STRATUM]["agreement_fraction"], 0)
        override = summarize_accuracy(rows, preparations=preps, tolerance_pp=0.625)
        self.assertEqual(override["tolerance_pp"], 0.625)
        self.assertEqual(override["targets_pp"], report["targets_pp"])
        self.assertTrue(all(summary["tolerance_pp"] == 0.625 for summary in override["summaries"].values()))
        self.assertTrue(all(summary["agreement_fraction"] == 1 for summary in override["summaries"].values()))

    def test_method_vendor_and_d1_strata(self):
        strata = [STRATUM, ("internal_standard", "JEOL", "long"),
                  ("external_calibration", "Bruker", "long"), ("external_calibration", "JEOL", "short")]
        preps = [preparation("a", stratum=key) for key in strata]
        rows = [result(prep, 99 + i, str(i)) for i, prep in enumerate(preps)]
        report = summarize_accuracy(rows, preparations=preps)
        self.assertEqual({key: item["bias_pp"] for key, item in report["summaries"].items()},
                         dict(zip(strata, range(4))))

    def test_abstentions_failures_missing_rows_and_predeclared_exclusions(self):
        a, b, c, d, e, f = (preparation(name) for name in "abcdef")
        rows = [result(a, 100), result(b, 99, status="fail"),
                result(c, None, status="cannot be assessed"), result(e, 90),
                result(f, 99, reference_purity_pct=None)]
        report = summarize_accuracy(rows, preparations=[a, b, c, d, e, f], exclusions={e: "pilot exposure"})
        summary = report["summaries"][STRATUM]
        self.assertEqual(summary["eligible_preparations"], 5)
        self.assertEqual(summary["result_preparations"], 1)
        self.assertEqual(summary["eligible_dependency_groups"], 4)
        self.assertEqual(summary["eligible_result_fraction"], 1 / 5)
        self.assertEqual(summary["bias_pp"], 1)
        self.assertEqual(report["preparations"][1]["statuses"]["fail"], 1)
        self.assertEqual(report["preparations"][2]["statuses"]["cannot be assessed"], 1)
        self.assertEqual(report["preparations"][4]["exclusion"], "pilot exposure")
        with self.assertRaisesRegex(ValueError, "declared preparation"):
            summarize_accuracy(rows, preparations=[a], exclusions={e: "pilot exposure"})
        with self.assertRaisesRegex(ValueError, "require.*a reason"):
            summarize_accuracy([], preparations=[a], exclusions={a: " "})

    def test_missing_evidence_markers(self):
        a = preparation("a")
        row = result(a, 100, reference_coverage_basis="not documented")
        report = summarize_accuracy([row], preparations=[a], k=2)
        prep = report["preparations"][0]
        summary = report["summaries"][STRATUM]
        self.assertIsNone(prep["standard_u_ratio"])
        self.assertIn("reference coverage basis", prep["error_u_reason"])
        self.assertEqual(summary["result_preparations"], 1)
        self.assertEqual(summary["standard_u_comparisons"], 0)
        self.assertIsNone(summary["standard_u_agreement_fraction"])
        self.assertIsNone(summary["bias_standard_u_pp"])
        with self.assertRaisesRegex(ValueError, "require.*a reason"):
            summarize_accuracy([row], preparations=[a], exclusions={a: "unknown"})

    def test_zero_denominators(self):
        a = preparation("a")
        for exclusions, result_fraction in (({}, 0), ({a: "pilot exposure"}, None)):
            with self.subTest(exclusions=exclusions):
                summary = summarize_accuracy([], preparations=[a], exclusions=exclusions)["summaries"][STRATUM]
                self.assertEqual(summary["eligible_result_fraction"], result_fraction)
                for key in ("bias_pp", "median_absolute_error_pp", "agreement_fraction", "bias_standard_u_pp",
                            "standard_u_agreement_fraction"):
                    self.assertIsNone(summary[key])
                self.assertEqual(summary["bias_u_reason"], "No eligible results")

    def test_normalization_requires_both_standard_uncertainties_basis_and_k(self):
        preps = [preparation(name) for name in "abcdef"]
        row = result(preps[0], 100, measurement_u_pp=0.3, reference_u_pp=0.4)
        rows = [row, result(preps[1], 99, reference_u_pp=None), result(preps[2], 99, measurement_u_pp=None),
                result(preps[3], 99, reference_coverage_basis="unknown"),
                result(preps[4], 99, measurement_u_pp=0, reference_u_pp=0),
                result(preps[5], 101, measurement_u_pp=0.3, reference_u_pp=0.4)]
        report = summarize_accuracy(rows, preparations=preps, k=2)
        self.assertEqual(report["uncertainty_comparison"], "standard-u comparison with factor k")
        self.assertEqual([entry["standard_u_ratio"] for entry in report["preparations"]],
                         [1, None, None, None, None, 2])
        summary = report["summaries"][STRATUM]
        self.assertEqual(summary["standard_u_agreement_fraction"], 0.5)
        self.assertEqual(summary["standard_u_comparisons"], 2)
        self.assertIsNone(summary["bias_standard_u_pp"])
        self.assertIn("standard uncertainty", summary["bias_u_reason"])
        report = summarize_accuracy([row], preparations=preps[:1])
        self.assertIsNone(report["preparations"][0]["standard_u_ratio"])

    def test_identical_rows_and_repeats_add_no_independent_units(self):
        a, b = preparation("a"), preparation("b")
        first, second = result(a, 99), result(b, 100)
        report = summarize_accuracy([first, second], preparations=[a, b])
        duplicate = summarize_accuracy([first, first, second], preparations=[a, b])
        self.assertEqual(report, duplicate)
        repeated = summarize_accuracy([first, replace(first, run_id="a2"), second], preparations=[a, b], k=2)
        self.assertEqual(repeated["preparations"][0]["eligible_repeats"], 2)
        self.assertEqual(repeated["preparations"][0]["repeat_sd_pp"], 0)
        self.assertIsNone(repeated["preparations"][0]["standard_u_ratio"])
        self.assertEqual(repeated["summaries"][STRATUM]["result_dependency_groups"], 2)
        with self.assertRaisesRegex(ValueError, "Conflicting rows"):
            summarize_accuracy([first, replace(first, status="fail")], preparations=[a])

    def test_repeats_can_use_different_calibrant_runs(self):
        a = preparation("a")
        first = result(a, 99, "a1", dependency_group="ec1", reference_u_pp=0)
        second = replace(first, run_id="a2", measured_purity_pct=101, dependency_group="ec2")
        failed = replace(first, run_id="a3", status="fail", dependency_group="ec3")
        rows = [first, second, failed]
        report = summarize_accuracy(rows, preparations=[a])
        summary = report["summaries"][STRATUM]
        self.assertEqual(len(report["preparations"]), 1)
        self.assertEqual(report["preparations"][0]["measured_purity_pct"], 100)
        self.assertEqual(summary["eligible_preparations"], 1)
        self.assertEqual(summary["result_preparations"], 1)
        self.assertEqual(summary["eligible_dependency_groups"], 3)
        self.assertEqual(summary["result_dependency_groups"], 2)
        self.assertAlmostEqual(summary["bias_standard_u_pp"], sqrt(0.5))
        with self.assertRaisesRegex(ValueError, "contradicts independent groups"):
            summarize_accuracy(rows, preparations=[a], error_covariances_pp2={frozenset(("a1", "a2")): 0.1})

    def test_shared_calibrant_propagates_weighted_error_covariance(self):
        a, b = preparation("a"), preparation("b")
        rows = [result(a, 99, "a1", dependency_group="calibrant"),
                result(a, 100, "a2", dependency_group="calibrant"),
                result(b, 98, dependency_group="calibrant")]
        # Each run has variance 0.75 residual + 0.25 reference + 0.25 calibrant.
        # Repeats share their reference, and all three runs share the calibrant.
        within = {frozenset(("a1", "a2")): 0.5}
        cross = {frozenset((run, "b")): 0.25 for run in ("a1", "a2")}
        shared = summarize_accuracy(rows, preparations=[a, b], error_covariances_pp2={**within, **cross})
        summary = shared["summaries"][STRATUM]
        self.assertEqual(summary["result_dependency_groups"], 1)
        self.assertEqual(summary["result_preparations"], 2)
        # Prep A variance is (1.25 + 1.25 + 2*0.5)/4.
        self.assertAlmostEqual(shared["preparations"][0]["error_standard_u_pp"], sqrt(0.875))
        # Bias variance is (0.875 + 1.25 + 2*0.25)/4 with equal prep weights.
        self.assertAlmostEqual(summary["bias_standard_u_pp"], sqrt(0.65625))
        independent = summarize_accuracy([*rows[:2], replace(rows[2], dependency_group="other calibrant")],
                                         preparations=[a, b], error_covariances_pp2=within)
        self.assertAlmostEqual(independent["summaries"][STRATUM]["bias_standard_u_pp"], sqrt(0.53125))
        self.assertEqual(independent["summaries"][STRATUM]["result_dependency_groups"], 2)
        missing = summarize_accuracy(rows, preparations=[a, b], error_covariances_pp2=within)["summaries"][STRATUM]
        self.assertIsNone(missing["bias_standard_u_pp"])
        self.assertIn("covariance", missing["bias_u_reason"])
        # Each pair is within its variance bound, but the full matrix is indefinite.
        inconsistent = {**dict.fromkeys(within, -1), **dict.fromkeys(cross, 1)}
        impossible = summarize_accuracy(rows, preparations=[a, b],
                                        error_covariances_pp2=inconsistent)["summaries"][STRATUM]
        self.assertIsNone(impossible["bias_standard_u_pp"])
        self.assertIn("positive semidefinite", impossible["bias_u_reason"])

    def test_zero_variance_requires_zero_covariance(self):
        a = preparation("a")
        rows = [result(a, 99, "a1", measurement_u_pp=0, reference_u_pp=0),
                result(a, 100, "a2", reference_u_pp=0)]
        pair = frozenset(("a1", "a2"))
        # These nonzero covariances are small enough to pass the matrix tolerance.
        for covariance in (1e-9, -1e-9):
            with self.subTest(covariance=covariance):
                report = summarize_accuracy(rows, preparations=[a], error_covariances_pp2={pair: covariance})
                summary = report["summaries"][STRATUM]
                self.assertIsNone(summary["bias_standard_u_pp"])
                self.assertIn("covariance exceeds the variance bound", summary["bias_u_reason"])
        report = summarize_accuracy(rows, preparations=[a], error_covariances_pp2={pair: 0})
        self.assertEqual(report["summaries"][STRATUM]["bias_standard_u_pp"], 0.5)

    def test_covariance_bound_with_positive_variances(self):
        a = preparation("a")
        rows = [result(a, 99, "a1", measurement_u_pp=0.7, reference_u_pp=0),
                result(a, 100, "a2", measurement_u_pp=0.1, reference_u_pp=0)]
        pair = frozenset(("a1", "a2"))
        at_bound = summarize_accuracy(rows, preparations=[a], error_covariances_pp2={pair: 0.07})
        self.assertAlmostEqual(at_bound["summaries"][STRATUM]["bias_standard_u_pp"], 0.4)
        above_bound = summarize_accuracy(rows, preparations=[a], error_covariances_pp2={pair: 0.0700000000001})
        self.assertIsNone(above_bound["summaries"][STRATUM]["bias_standard_u_pp"])

    def test_perfect_correlation_survives_rounding_of_the_bound(self):
        # sqrt(0.01**2 * 0.37**2) rounds below 0.01 * 0.37; the product itself must pass.
        a = preparation("a")
        rows = [result(a, 99, "a1", measurement_u_pp=0.01, reference_u_pp=0),
                result(a, 100, "a2", measurement_u_pp=0.37, reference_u_pp=0)]
        report = summarize_accuracy(rows, preparations=[a],
                                    error_covariances_pp2={frozenset(("a1", "a2")): 0.01 * 0.37})
        summary = report["summaries"][STRATUM]
        self.assertIsNone(summary["bias_u_reason"])
        self.assertAlmostEqual(summary["bias_standard_u_pp"], (0.01 + 0.37) / 2)


class ReleaseGateTests(unittest.TestCase):
    def evaluate(self, rows, preps, **changes):
        options = {"k": 2, "k_basis": "Documented coverage factor", **changes}
        return evaluate_release_gates(rows, preparations=preps, **options)[STRATUM[0]]

    def test_all_gates_pass_with_seven_of_eight_below_target(self):
        preps = [preparation(str(i)) for i in range(8)]
        rows = [result(prep, 99.75 if i == 7 else 99) for i, prep in enumerate(preps)]
        report = self.evaluate(rows, preps)
        gates = report["gates"]
        self.assertEqual(report["status"], "pass")
        self.assertEqual({gate["status"] for gate in gates.values()}, {"pass"})
        sign = gates["sign_test"]
        self.assertEqual((sign["successes"], sign["failures"], sign["units"]), (7, 1, 8))
        self.assertAlmostEqual(sign["p_value"], 0.03515625)
        self.assertEqual(gates["accuracy"]["median_absolute_error_pp"], 0)

    def test_method_targets_apply_independently_at_the_accuracy_boundary(self):
        preps = [preparation(str(i), key) for key in (STRATUM, INTERNAL_STRATUM) for i in range(8)]
        rows = [result(prep, 99.75, str(i)) for i, prep in enumerate(preps)]
        report = evaluate_release_gates(rows, preparations=preps, k=2, k_basis="Documented coverage factor")
        external = report[STRATUM[0]]["gates"]
        internal = report[INTERNAL_STRATUM[0]]["gates"]
        self.assertEqual(external["accuracy"]["target_pp"], 0.75)
        self.assertEqual(external["accuracy"]["status"], "pass")
        self.assertEqual(external["sign_test"]["status"], "fail")
        self.assertEqual(internal["accuracy"]["target_pp"], 0.50)
        self.assertEqual(internal["accuracy"]["status"], "fail")
        self.assertEqual(internal["independent_units"]["scored_units"], 8)
        self.assertEqual(report[INTERNAL_STRATUM[0]]["status"], "fail")

    def test_shared_group_chain_spans_strata_and_exclusions_cannot_bridge(self):
        a, b, c, d, e, f, g, excluded = [preparation(name) for name in "abcdefgh"]
        b = replace(b, vendor="Bruker")
        c = replace(c, d1_stratum="short")
        sibling = replace(excluded, d1_stratum="short")
        rows = [result(a, 99.25, dependency_group="left"),
                result(b, 98, "b1", dependency_group="left"), result(b, 100, "b2", dependency_group="right"),
                result(c, 99.25, dependency_group="right"), *[result(prep, 99) for prep in (d, e, f, g)],
                result(sibling, 99, "h1", dependency_group="right"), result(sibling, 99, "h2", dependency_group="d")]
        report = self.evaluate(rows, [a, b, c, d, e, f, g, excluded, sibling], exclusions={excluded: "pilot exposure"})
        gates = report["gates"]
        self.assertEqual(report["status"], "open")
        self.assertEqual(gates["accuracy"]["median_absolute_error_pp"], 0)
        independence = gates["independent_units"]
        self.assertEqual((independence["status"], independence["scored_units"], independence["eligible_units"]),
                         ("open", 5, 5))
        self.assertEqual(gates["sign_test"]["successes"], 5)
        self.assertEqual(gates["sign_test"]["p_value"], 1 / 32)
        self.assertEqual(gates["fraction_assessed"]["eligible_preparations"], 7)

    def test_delay_strata_are_pooled_before_scoring_physical_preparations(self):
        long = [preparation(str(i)) for i in range(8)]
        short = [replace(prep, d1_stratum="short") for prep in long]
        rows = [result(prep, 98, f"long-{i}") for i, prep in enumerate(long)]
        rows.extend(result(prep, 100.5, f"short-{i}", dependency_group=f"short-{i}") for i, prep in enumerate(short))
        report = self.evaluate(rows, [*long, *short, preparation("no runs")])
        gates = report["gates"]
        assessed = gates["fraction_assessed"]
        self.assertEqual((report["status"], assessed["status"], assessed["fraction"]), ("fail", "fail", 8 / 9))
        self.assertEqual(gates["accuracy"]["median_absolute_error_pp"], 0.25)
        self.assertEqual(gates["independent_units"]["scored_units"], 8)
        self.assertEqual((gates["sign_test"]["successes"], gates["sign_test"]["units"]), (8, 9))
        self.assertEqual(gates["uncertainty_agreement"]["compared_preparations"], 8)

    def test_passing_repeats_do_not_mask_flags_or_abstentions_in_sign_test(self):
        preps = [preparation(str(i)) for i in range(8)]
        passing = [result(prep, 99) for prep in preps]
        for status in ("fail", "cannot be assessed"):
            with self.subTest(status=status):
                repeats = [result(prep, None, f"repeat-{i}", status=status) for i, prep in enumerate(preps)]
                gates = self.evaluate([*passing, *repeats], preps)["gates"]
                sign = gates["sign_test"]
                self.assertEqual((sign["status"], sign["successes"], sign["units"]), ("fail", 0, 8))
                self.assertEqual(gates["fraction_assessed"]["scored_preparations"], 8)

    def test_passing_repeats_with_missing_purities_are_sign_failures(self):
        preps = [preparation(str(i)) for i in range(8)]
        passing = [result(prep, 99) for prep in preps]
        baseline = self.evaluate(passing, preps)["gates"]
        for field in ("measured_purity_pct", "reference_purity_pct"):
            with self.subTest(field=field):
                repeats = [replace(row, run_id=f"repeat-{row.run_id}", **{field: None}) for row in passing]
                gates = self.evaluate([*passing, *repeats], preps)["gates"]
                sign = gates["sign_test"]
                self.assertEqual((sign["status"], sign["successes"], sign["failures"], sign["units"]), ("fail", 0, 8, 8))
                self.assertEqual(gates["accuracy"], baseline["accuracy"])
                self.assertEqual(gates["fraction_assessed"], baseline["fraction_assessed"])

    def test_unscored_preparations_are_sign_failures_and_stay_in_denominators(self):
        preps = [preparation(str(i)) for i in range(12)]
        rows = [result(prep, 99) for prep in preps[:8]]
        rows.extend([result(preps[8], 99, status="fail"), result(preps[9], None, status="cannot be assessed"),
                     result(preps[10], 99, reference_purity_pct=None)])
        gates = self.evaluate(rows, preps)["gates"]
        self.assertEqual(gates["independent_units"]["scored_units"], 8)
        sign = gates["sign_test"]
        self.assertEqual((sign["status"], sign["successes"], sign["failures"], sign["units"]), ("fail", 8, 4, 12))
        self.assertEqual(sign["p_value"], 794 / 4096)
        assessed = gates["fraction_assessed"]
        self.assertEqual((assessed["status"], assessed["fraction"]), ("fail", 8 / 12))

    def test_fraction_assessed_accepts_nine_of_ten(self):
        preps = [preparation(str(i)) for i in range(10)]
        rows = [result(prep, 99) for prep in preps[:9]]
        rows.append(result(preps[9], 99, status="fail", dependency_group="0"))
        report = self.evaluate(rows, preps)
        self.assertEqual(report["status"], "pass")
        gates = report["gates"]
        self.assertEqual(gates["fraction_assessed"]["fraction"], 0.90)
        self.assertEqual(gates["independent_units"]["scored_units"], 9)
        self.assertEqual(gates["sign_test"]["successes"], 8)
        self.assertEqual(gates["sign_test"]["units"], 9)

    def test_uncertainty_accepts_boundary_and_lists_every_miss(self):
        preps = [preparation(str(i)) for i in range(10)]
        rows = [result(prep, 99, measurement_u_pp=0, reference_u_pp=0) for prep in preps[:8]]
        rows.extend([result(preps[8], 99.5, measurement_u_pp=0.25, reference_u_pp=0),
                     result(preps[9], 99.625, measurement_u_pp=0.25, reference_u_pp=0)])
        gate = self.evaluate(rows, preps)["gates"]["uncertainty_agreement"]
        self.assertEqual((gate["status"], gate["fraction"], gate["misses"]), ("pass", 0.90, ["9"]))
        rows[8] = replace(rows[8], measurement_u_pp=0)
        report = self.evaluate(rows, preps)
        gate = report["gates"]["uncertainty_agreement"]
        self.assertEqual((gate["status"], gate["fraction"], gate["misses"]), ("fail", 0.80, ["8", "9"]))
        self.assertEqual(report["status"], "fail")

    def test_uncertainty_opens_with_reasons_and_keeps_observed_misses(self):
        preps = [preparation(str(i)) for i in range(8)]
        rows = [result(prep, 99) for prep in preps]
        for options, reason in (({"k": None}, "factor k"), ({"k_basis": None}, "basis k_basis"),
                                ({"k_basis": "not documented"}, "basis k_basis")):
            with self.subTest(options=options):
                report = self.evaluate(rows, preps, **options)
                gate = report["gates"]["uncertainty_agreement"]
                self.assertEqual((report["status"], gate["status"]), ("open", "open"))
                self.assertIn(reason, gate["reason"])
        rows[0] = replace(rows[0], reference_coverage_basis=None)
        rows[1] = replace(rows[1], measured_purity_pct=99.625, measurement_u_pp=0, reference_u_pp=0)
        rows.append(replace(rows[2], run_id="2-repeat"))
        report = self.evaluate(rows, preps)
        gate = report["gates"]["uncertainty_agreement"]
        self.assertEqual((report["status"], gate["status"], gate["misses"]), ("open", "open", ["1"]))
        self.assertEqual(gate["compared_preparations"], 6)
        self.assertIsNone(gate["fraction"])
        missing = {entry["preparation_id"]: entry["reason"] for entry in gate["missing_uncertainty"]}
        self.assertEqual(set(missing), {"0", "2"})
        self.assertIn("reference coverage basis", missing["0"])
        self.assertIn("covariance", missing["2"])

    def test_empty_eligible_roster_leaves_all_gates_open(self):
        excluded = preparation("pilot")
        report = self.evaluate([result(excluded, 90)], [excluded], exclusions={excluded: "pilot exposure"})
        self.assertEqual(report["status"], "open")
        self.assertEqual({gate["status"] for gate in report["gates"].values()}, {"open"})
        self.assertTrue(all(gate["reason"] for gate in report["gates"].values()))


if __name__ == "__main__":
    unittest.main()
