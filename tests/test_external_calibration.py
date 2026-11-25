import json
import unittest
from dataclasses import replace
from fractions import Fraction
from math import isclose, nan
from os import scandir
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from spinivo.acquisition import RecoveryEvidence, SignalEvidence, check_recovery
from spinivo.external_calibration import (
    ExternalCalibrationAssay, ExternalRun, PairRecord, _calculate_external_calibration,
    calculate_external_calibration, link_external_checks,
)
from spinivo.external_calibration_report import purity_report
from spinivo.integral_components import ComponentEvidence, Line, check_region
from spinivo.internal_standard import CALCULATION_ONLY, Correlation, CoverageFactor, Prerequisite, Quantity
from spinivo.processing import QUANTITATIVE_SETTINGS, process
from spinivo.readers.bruker import read_bruker
from spinivo.readers.jeol import read_jeol
from tests.integral_fixtures import NOISE_BOUNDS_HZ, ppm_bounds
from tests.test_internal_standard import NOISE_GAP, acquisition_checks, clean_regions


PREPARATION = ("m_s", "w_s", "M_s", "V_s")
QUANTITIES = (*PREPARATION, "M_a", "V_a", "m_a", "K", "D")
RECORDED = ("relaxation_delay_s", "acquisition_time_s", "spectral_width_hz", "complex_points", "observe_frequency_mhz")


def load_assay_case(transfer=True):
    """Return the FID02 pair, by default with the fixture's illustrative transfer evidence."""
    case = json.loads((Path(__file__).parent / "fixtures" / "external_calibration_cases.json").read_text())
    inputs = case["assay"]
    for name in ("sample", "calibrant"):
        inputs[name] = ExternalRun(**inputs[name])
    for name in QUANTITIES:
        fields = inputs[name]
        inputs[name] = Quantity(**{**fields, "source_ids": tuple(fields.get("source_ids", ()))})
    for name in ("sample_signals", "calibrant_signals", "R_source_ids"):
        inputs[name] = tuple(inputs.get(name, ()))
    inputs["correlations"] = tuple(Correlation(tuple(c["inputs"]), c["r"], c["source"]) for c in inputs["correlations"])
    inputs["k"] = CoverageFactor(**inputs["k"])
    assay = ExternalCalibrationAssay(**inputs)
    if not transfer:
        return assay, case["expected"]
    extra = case["illustrative_transfer"]
    runs = {role: replace(getattr(assay, role), **extra["runs"]) for role in ("sample", "calibrant")}
    assay = replace(assay, pair=PairRecord(**extra["pair"]), S=Quantity(**extra["S"]), S_gains=tuple(extra["S_gains"]),
                    S_storage_scales=tuple(extra["S_storage_scales"]), **runs)
    return assay, {**case["expected"], **extra["expected"]}


class ExternalCalibrationTests(unittest.TestCase):
    def setUp(self):
        self.assay, self.expected = load_assay_case()
        self.original = calculate_external_calibration(self.assay)
        self.areas = (self.original.J_a, self.original.J_s, self.original.R)

    def assertRelative(self, actual, expected):
        self.assertTrue(isclose(actual, expected, rel_tol=1e-9), (actual, expected))

    def assertBlocked(self, assay, reason, areas):
        result = calculate_external_calibration(assay)
        self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
        self.assertIn(reason, result.reasons)
        self.assertEqual((result.J_a, result.J_s, result.R), areas)

    def assertTransferBlocked(self, assay, reason):
        self.assertBlocked(assay, reason, self.areas)
        self.assertEqual(calculate_external_calibration(assay).c_s, self.original.c_s)

    def test_hand_worked_budget(self):
        result, expected = self.original, self.expected
        self.assertEqual((result.status, result.reasons, result.concentration_unit),
                         ("cannot be assessed", (CALCULATION_ONLY,), "mmol/L"))
        self.assertEqual(self.areas[:2], (0.55648986 / 32, 0.35070252 / 32))
        for name in ("R", "c_s", "c_a"):
            self.assertRelative(getattr(result, name), expected[name])
        self.assertEqual(self.assay.K.value, 7.36572 / 7.37538)
        self.assertRelative(self.assay.K.value, expected["K"])
        self.assertRelative(result.value, expected["purity"])
        self.assertRelative(result.value / 0.999, expected["purity_relative_to_sheet_basis"])
        budget = result.budget
        self.assertEqual(budget.k, self.assay.k)
        for name in ("variance", "standard_u", "relative_standard_u", "expanded_u"):
            self.assertRelative(getattr(budget, name), expected[name])
        self.assertEqual({item.name for item in budget.inputs}, set(expected["sensitivities"]))
        for item in budget.inputs:
            self.assertIsInstance(item.quantity.source_ids, tuple)
            self.assertRelative(item.sensitivity, expected["sensitivities"][item.name])
        self.assertEqual(tuple(item.correlation for item in budget.covariances), self.assay.correlations)
        for item, term in zip(budget.covariances, expected["covariance_terms"]):
            self.assertRelative(item.term, term)

    def test_direct_calculation_without_run_preserves_independent_values(self):
        for role, areas in (("sample", (None, self.original.J_s)), ("calibrant", (self.original.J_a, None)),
                            ("both", (None, None))):
            with self.subTest(role=role):
                missing = ("sample", "calibrant") if role == "both" else (role,)
                assay = replace(self.assay, **dict.fromkeys(missing))
                result = calculate_external_calibration(assay)
                self.assertEqual((result.status, result.value, result.budget, result.R, result.c_a),
                                 ("cannot be assessed", None, None, None, None))
                self.assertEqual((result.J_a, result.J_s, result.c_s), (*areas, self.original.c_s))
                report = purity_report(result)
                for name in missing:
                    self.assertIn(f"{name} run: Missing run record", report["reasons"])

    def test_non_unity_S_scales_the_sheet_budget(self):
        _, sheet = load_assay_case(transfer=False)
        S, u, expected = self.assay.S.value, self.assay.S.standard_u, self.expected
        self.assertEqual(S, 0.98)
        for name in ("c_a", "purity", "purity_relative_to_sheet_basis", "purity_at_K_1", "purity_at_K_1_02"):
            self.assertRelative(expected[name], S * sheet[name])
        self.assertRelative(expected["variance"], S**2 * sheet["variance"] + (sheet["purity"] * u)**2)
        for name, value in sheet["sensitivities"].items():
            self.assertRelative(expected["sensitivities"][name], S * value)
        self.assertRelative(expected["sensitivities"]["S"], sheet["purity"])
        for term, value in zip(expected["covariance_terms"], sheet["covariance_terms"]):
            self.assertRelative(term, S**2 * value)
        unity = calculate_external_calibration(replace(self.assay, S=replace(self.assay.S, value=1)))
        self.assertRelative(unity.value, sheet["purity"])

    def test_scan_normalisation_and_corrections(self):
        a, s = self.assay.sample, self.assay.calibrant
        doubled = replace(self.assay, sample=replace(a, area=2 * a.area, scans=2 * a.scans))
        averaged = replace(doubled, calibrant=replace(s, area=s.area / s.scans, area_mode="averaged"))
        both = replace(averaged, sample=replace(a, area=a.area / a.scans, area_mode="averaged"))
        reader = replace(self.assay, calibrant=replace(s, area=4 * s.area, reader_scale=4))
        for assay in (doubled, averaged, both, reader):
            self.assertEqual(calculate_external_calibration(assay), self.original)
        self.assertBlocked(replace(both, sample=replace(both.sample, scans=None)), "sample scans is missing", self.areas)
        unity, adjusted = (calculate_external_calibration(replace(self.assay, K=replace(self.assay.K, value=value)))
                           for value in (1, 1.02))
        self.assertRelative(unity.value, self.expected["purity_at_K_1"])
        self.assertRelative(adjusted.value, self.expected["purity_at_K_1_02"])
        self.assertEqual(adjusted.value / unity.value, 1.02)
        drifted = calculate_external_calibration(replace(self.assay, D=replace(self.assay.D, value=0.98)))
        self.assertRelative(drifted.value / self.original.value, 0.98)

    def test_certified_concentration(self):
        stock = ("stock:fixture-certified",)
        run = ("calibrant-run:" + self.assay.calibrant.run_id,)
        c_s = Quantity(self.original.c_s, 0.01, "mmol/L", "Illustrative certificate", stock)
        zero = Correlation(("R", "D"), 0, "Illustrative drift validation", "Separate residual estimates")
        certified = replace(self.assay, c_s_route="certified", c_s=c_s, **dict.fromkeys(PREPARATION),
                            R_source_ids=run, D=replace(self.assay.D, source_ids=stock + run),
                            correlations=(Correlation(("c_s", "D"), 0.2, "Illustrative correlation"), zero))
        result = calculate_external_calibration(certified)
        self.assertEqual((result.value, result.c_a), (self.original.value, self.original.c_a))
        self.assertEqual(result.budget.covariances[-1].term, 0)
        self.assertBlocked(replace(certified, correlations=certified.correlations[1:]),
                           "missing covariance for c_s, D (shared source stock:fixture-certified)", self.areas)
        self.assertBlocked(replace(certified, correlations=certified.correlations[:1]),
                           f"missing covariance for R, D (shared source {run[0]})", self.areas)
        self.assertEqual({item.name for item in result.budget.inputs},
                         {"R", "c_s", "M_a", "V_a", "m_a", "K", "D", "S"})
        for name, wrong, reason in (
            ("c_s", None, "c_s is missing"),
            ("c_s", replace(c_s, value=None), "c_s.value is missing"),
            ("c_s", replace(c_s, standard_u=None), "c_s.standard_u is missing"),
            ("c_s", replace(c_s, unit="mol/L"), "c_s.unit must be mmol/L"),
            ("m_a", replace(certified.m_a, unit="g"), "m_a.unit must be mg for certified c_s"),
            ("V_a", replace(certified.V_a, unit="L"), "V_a.unit must be mL for certified c_s"),
            ("M_a", replace(certified.M_a, unit="kg/mol"), "M_a.unit must be g/mol for certified c_s"),
        ):
            with self.subTest(reason=reason):
                invalid = replace(certified, **{name: wrong})
                self.assertBlocked(invalid, reason, self.areas)
                self.assertEqual(calculate_external_calibration(replace(invalid, **{name: getattr(certified, name)})), result)
        wrong_route = replace(self.assay, c_s_route="certified")
        self.assertBlocked(wrong_route, "c_s is missing", self.areas)
        self.assertIsNone(calculate_external_calibration(wrong_route).c_s)

    def test_run_blockers_and_restoration(self):
        conditions = ("instrument", "probe", "pulse_sequence", "decoupling", "excitation", "temperature", "solvent",
                      "tube", "fill")
        for role, incomplete in (("sample", (None, self.areas[1], None)), ("calibrant", (self.areas[0], None, None))):
            run = getattr(self.assay, role)
            cases = [(name, None, f"{name} is missing", self.areas)
                     for name in ("run_id", "preparation_id", "receiver_gain", *conditions)]
            cases += [(name, None, f"{name} is missing", incomplete) for name in ("scans", "area_mode", "reader_scale")]
            cases += [("scans", value, "scans must be finite and positive", incomplete) for value in (0, -1, True)]
            cases += [("scans", 1.5, "scans must be a whole count", incomplete),
                      ("reader_scale", 0, "reader_scale must be finite and positive", incomplete),
                      ("area_mode", 42, "area_mode must be summed or averaged", incomplete)]
            for name, value, reason, areas in cases:
                with self.subTest(role=role, field=name, value=value):
                    changed = replace(self.assay, **{role: replace(run, **{name: value})})
                    self.assertBlocked(changed, f"{role} {reason}", areas)
                    self.assertEqual(calculate_external_calibration(replace(changed, **{role: run})), self.original)
        values = ("other", "other probe", "other sequence", "no decoupling", "30 degree pulse", 26, "water", "3 mm",
                  "half fill")
        for name, value in zip(conditions, values):
            with self.subTest(mismatch=name):
                changed = replace(self.assay, sample=replace(self.assay.sample, **{name: value}))
                if name == "instrument":
                    self.assertTransferBlocked(changed,
                                               "instrument differs; transfer between instruments is not supported")
                    self.assertTransferBlocked(replace(changed, K_covers=(name,)), "K cannot cover instrument")
                    continue
                self.assertTransferBlocked(changed, f"{name} differs without K evidence")
                self.assertBlocked(replace(changed, K_covers=name), f"{name} differs without K evidence", self.areas)
                covered = replace(changed, K_covers=(name,))
                self.assertEqual(calculate_external_calibration(covered), self.original)
                missing = replace(covered, sample=replace(covered.sample, **{name: None}))
                self.assertBlocked(missing, f"sample {name} is missing", self.areas)
                self.assertEqual(calculate_external_calibration(replace(changed, sample=self.assay.sample)), self.original)

    def test_tolerances_and_recorded_settings(self):
        recorded = dict(zip(RECORDED, (20.0, 8.745124, 7494.005, 65536, 399.782198)))
        base = replace(self.assay, **{role: replace(getattr(self.assay, role), **recorded)
                                      for role in ("sample", "calibrant")})
        for changes, reason in (
            ({"temperature": 25.05}, None), ({"temperature": 24.8}, "temperature differs without K evidence"),
            ({"relaxation_delay_s": 60.0, "complex_points": 32768, "spectral_width_hz": 6000.0}, None),
            ({"observe_frequency_mhz": 399.79}, None),
            ({"observe_frequency_mhz": 500.13},
             "observe_frequency_mhz differs; transfer between instruments is not supported"),
        ):
            with self.subTest(changes=changes):
                changed = replace(base, sample=replace(base.sample, **changes))
                if reason is None:
                    self.assertEqual(calculate_external_calibration(changed), self.original)
                else:
                    self.assertTransferBlocked(changed, reason)

    def test_pair_record_is_required(self):
        pair = self.assay.pair
        reused = "pair source names a run file or folder, not a pair record"
        cases = [(None, "pair record is missing"), ((), "pair must be a PairRecord"),
                 (replace(pair, kind="guess"), "pair kind must be record or attestation"),
                 (replace(pair, source=" "), "pair source is missing"),
                 (replace(pair, source=self.assay.sample.run_id), reused), (replace(pair, source="."), reused)]
        for role in ("sample", "calibrant"):
            for name in ("run_id", "preparation_id"):
                field, record = f"{role}_{name}", getattr(getattr(self.assay, role), name)
                differs = f"pair {field} 'x' differs from the {role} record {record!r}"
                cases += [(replace(pair, **{field: None}), f"pair {field} is missing"),
                          (replace(pair, **{field: "x"}), differs)]
        for value, reason in cases:
            with self.subTest(reason=reason):
                changed = replace(self.assay, pair=value)
                self.assertTransferBlocked(changed, reason)
                self.assertEqual(calculate_external_calibration(replace(changed, pair=pair)), self.original)
        attested = replace(self.assay, pair=replace(pair, kind="attestation", source="Analyst attestation, signed log"))
        self.assertEqual(calculate_external_calibration(attested), self.original)

    def test_scale_evidence_is_required(self):
        S, K = self.assay.S, self.assay.K
        unequal = replace(self.assay, sample=replace(self.assay.sample, receiver_gain=31))
        cases = (
            (replace(self.assay, S=None), ("S is missing", "equal nominal receiver gains are not scale evidence")),
            (replace(unequal, S=None), ("S is missing", "receiver gains differ without a sourced S correction")),
            (unequal, ("S_gains (30, 30) do not cover the run receiver gains (31, 30)",)),
            (replace(unequal, K_covers=("receiver_gain",)), ("K cannot cover receiver_gain",)),
            (replace(self.assay, S_gains=None), ("S_gains is missing",)),
            (replace(self.assay, S=replace(S, value=1, standard_u=0)), ("S.standard_u must be finite and positive",)),
            (replace(self.assay, S=replace(S, source=" ")), ("S.source is missing",)),
            (replace(self.assay, S=replace(S, unit="mg")), ("S.unit must be dimensionless",)),
            (replace(self.assay, S=replace(S, source_ids=("run:x",)), K=replace(K, source_ids=("run:x",))),
             ("missing covariance for K, S (shared source run:x)",)),
        )
        for assay, reasons in cases:
            for reason in reasons:
                with self.subTest(reason=reason):
                    self.assertTransferBlocked(assay, reason)
        self.assertEqual(calculate_external_calibration(replace(unequal, S_gains=(31, 30))), self.original)

    def test_pair_source_cannot_be_an_absolute_run_path_or_parent(self):
        run = Path(self.assay.sample.run_id)
        for source in (Path.cwd(), run.absolute(), Path.cwd() / "docs" / ".." / run):
            with self.subTest(source=source):
                assay = replace(self.assay, pair=replace(self.assay.pair, source=str(source)))
                self.assertTransferBlocked(assay, "pair source names a run file or folder, not a pair record")
        for folder in (Path("spectra/sample/10"), Path(__file__).parent / "fixtures/bruker/zenodo17083860_833-1/10"):
            moved = replace(self.assay, sample=replace(self.assay.sample, run_id=str(folder)),
                            pair=replace(self.assay.pair, sample_run_id=str(folder)))
            for name in ("fid", "acqus", "acqu", "pulseprogram", "pdata/1/procs", "../10/pdata", "../pair record"):
                with self.subTest(folder=folder, name=name):
                    assay = replace(moved, pair=replace(moved.pair, source=str(folder / name)))
                    if name == "../pair record":
                        self.assertEqual(calculate_external_calibration(assay), self.original)
                    else:
                        self.assertTransferBlocked(assay, "pair source names a run file or folder, not a pair record")

    def test_pair_source_must_be_independent_of_both_run_folders(self):
        root = Path(self.enterContext(TemporaryDirectory()))
        unrelated = root / "unrelated"
        unrelated.mkdir()
        for role in ("sample", "calibrant"):
            folder = root / role
            deep = folder / "pdata/1/a/b/entry.txt"
            deep.parent.mkdir(parents=True)
            deep.touch()
            alias = unrelated / f"{role}.txt"
            alias.hardlink_to(deep)
            for source in (unrelated, unrelated / "other.JDF", unrelated / "fid", unrelated / "ser", deep, alias):
                with self.subTest(source=source):
                    self.assertTransferBlocked(self.pair_assay(folder, source, role),
                                               "pair source names a run file or folder, not a pair record")

    def pair_assay(self, folder, source, role="sample"):
        run = replace(getattr(self.assay, role), run_id=str(folder))
        pair = replace(self.assay.pair, source=str(source), **{f"{role}_run_id": str(folder)})
        return replace(self.assay, pair=pair, **{role: run})

    def test_unreadable_run_subdirectory_blocks_pair_evidence(self):
        with TemporaryDirectory() as root:
            hidden, source = Path(root) / "run/private", Path(root) / "pair.txt"
            hidden.mkdir(parents=True)
            source.touch()
            (hidden / "entry.txt").hardlink_to(source)

            def scan(path):
                if Path(path) == hidden:
                    raise PermissionError(13, "Permission denied", str(hidden))
                return scandir(path)

            with patch("os.scandir", side_effect=scan):
                reason = f"pair source or run folder could not be inspected: [Errno 13] Permission denied: {str(hidden)!r}"
                self.assertTransferBlocked(self.pair_assay(hidden.parent, source), reason)

    def test_pair_record_can_be_beside_a_run_folder_named_fid_or_ser(self):
        with TemporaryDirectory() as root:
            source = Path(root) / "pair.txt"
            source.touch()
            for name in ("fid", "ser"):
                folder = Path(root) / name
                folder.mkdir()
                self.assertEqual(calculate_external_calibration(self.pair_assay(folder, source)), self.original)

    def test_each_run_needs_valid_acquisition_settings(self):
        for role in ("sample", "calibrant"):
            for field, value, reason in (
                ("observe_frequency_mhz", nan, "must be finite and positive"),
                ("observe_frequency_mhz", None, "is missing"),
                ("relaxation_delay_s", None, "is missing"),
                ("acquisition_time_s", -1, "must be finite and positive"),
                ("spectral_width_hz", 0, "must be finite and positive"),
                ("complex_points", 1.5, "must be a whole count"),
            ):
                with self.subTest(role=role, field=field, value=value):
                    run = replace(getattr(self.assay, role), **{field: value})
                    self.assertTransferBlocked(replace(self.assay, **{role: run}), f"{role} {field} {reason}")

    def test_placeholders_cannot_supply_matching_transfer_evidence(self):
        first, second = self.assay.correlations
        for marker in (" UNKNOWN ", "NOT FOUND", "n/a", "not_recorded", "NF"):
            with self.subTest(marker=marker):
                runs = {role: replace(getattr(self.assay, role), probe=marker) for role in ("sample", "calibrant")}
                assay = replace(self.assay, pair=replace(self.assay.pair, source=marker),
                                S=replace(self.assay.S, source=marker), correlations=(replace(first, source=marker), second), **runs)
                for reason in ("sample probe is missing", "calibrant probe is missing",
                               "pair source is missing", "S.source is missing", "correlation m_s, m_a source is missing"):
                    self.assertTransferBlocked(assay, reason)

    def test_temperature_limit_and_unrepresentable_reader_scale(self):
        for temperature, reason in ((25.1, None), (25.100001, "temperature differs without K evidence")):
            with self.subTest(temperature=temperature):
                assay = replace(self.assay, sample=replace(self.assay.sample, temperature=temperature))
                if reason:
                    self.assertTransferBlocked(assay, reason)
                else:
                    self.assertEqual(calculate_external_calibration(assay), self.original)
        sample = replace(self.assay.sample, reader_scale=Fraction(1, 10**400))
        self.assertBlocked(replace(self.assay, sample=sample), "sample reader_scale is outside the numeric range",
                           (None, self.areas[1], None))

    def test_fixture_pair_without_transfer_evidence_is_blocked(self):
        assay, _ = load_assay_case(transfer=False)
        result = calculate_external_calibration(assay)
        self.assertEqual((result.status, result.value, result.c_s), ("cannot be assessed", None, self.original.c_s))
        self.assertEqual((result.J_a, result.J_s, result.R), (None, None, None))
        missing = [f"{role} {name} is missing" for role in ("sample", "calibrant")
                   for name in ("instrument", "excitation", "fill", "reader_scale")]
        for reason in ("pair record is missing", "S is missing", "equal nominal receiver gains are not scale evidence",
                       *missing):
            self.assertIn(reason, result.reasons)

    def test_assay_blockers_and_restoration(self):
        first, second = self.assay.correlations
        cases = [(name, None, f"{name} is missing") for name in QUANTITIES]
        cases.append(("c_s_route", None, "c_s_route must be declared as certified or prepared"))
        cases += [(name, replace(getattr(self.assay, name), standard_u=None), f"{name}.standard_u is missing")
                  for name in PREPARATION]
        cases.append(("K_covers", ("pulse_sequnce",), "K_covers contains unknown condition 'pulse_sequnce'"))
        cases += [(name, replace(getattr(self.assay, name), source=""), f"{name}.source is missing") for name in ("K", "D")]
        cases += [(name, replace(getattr(self.assay, name), unit="other"), f"{name} and {sample} must use the same unit")
                  for name, sample in (("m_s", "m_a"), ("V_s", "V_a"), ("M_s", "M_a"))]
        cases.append(("c_s", Quantity(2.7, 0.01, "mmol/L", "Illustrative certificate"),
                      "only the declared prepared c_s inputs may be supplied"))
        cases.extend(("prerequisites", (Prerequisite("sample integral", status, "overlap"),), "sample integral: overlap")
                     for status in ("flag", "cannot be assessed"))
        cases.extend(("correlations", correlations, reason) for correlations, reason in (
            ((None,), "correlations must contain Correlation records or Covariance records"),
            ((replace(first, r=1.1),), "correlation m_s, m_a requires a finite r in [-1, 1]"),
            ((replace(second, source=""),), "correlation V_s, V_a source is missing"),
            ((replace(first, inputs=("R", "c_s")),), "correlation R, c_s contains an unknown budget input"),
            ((first, replace(first, inputs=("m_a", "m_s"))), "correlation m_a, m_s is a duplicate pair"),
            ((second,), "missing covariance for m_s, m_a (shared source balance:fixture-FID02)"),
            ((first,), "missing covariance for V_s, V_a (shared source volume:fixture-FID02)"),
            ((replace(first, r=1), second, Correlation(("R", "m_s"), 1, "Illustrative")),
             "correlation matrix is not positive semidefinite for R, m_s, m_a"),
        ))
        for name, value, reason in cases:
            with self.subTest(input=name, reason=reason):
                changed = replace(self.assay, **{name: value})
                self.assertBlocked(changed, reason, self.areas)
                self.assertEqual(calculate_external_calibration(replace(changed, **{name: getattr(self.assay, name)})),
                                 self.original)

    def test_identical_run_ids_block(self):
        calibrant = replace(self.assay.calibrant, run_id=self.assay.sample.run_id)
        self.assertBlocked(replace(self.assay, calibrant=calibrant), "sample and calibrant run_id must differ", self.areas)

    def test_nonpositive_areas_keep_signed_per_scan_values(self):
        a, s = self.assay.sample, self.assay.calibrant
        for mode in ("summed", "averaged"):
            for J_a, J_s, R, role in ((0, 2, 0, "sample"), (-1, 2, -0.5, "sample"), (1, -2, -0.5, "calibrant"),
                                      (-1, -2, 0.5, "sample"), (1, 0, None, "calibrant")):
                with self.subTest(mode=mode, J_a=J_a, J_s=J_s):
                    sample = replace(a, area=J_a * a.scans if mode == "summed" else J_a, area_mode=mode)
                    calibrant = replace(s, area=J_s * s.scans if mode == "summed" else J_s, area_mode=mode)
                    self.assertBlocked(replace(self.assay, sample=sample, calibrant=calibrant),
                                       f"{role} area must be finite and positive", (J_a, J_s, R))

    def test_nonzero_areas_outside_numeric_range_block(self):
        for mode, area in (("averaged", Fraction(1, 10**400)), ("summed", 5e-324)):
            with self.subTest(mode=mode):
                calibrant = replace(self.assay.calibrant, area=area, area_mode=mode)
                self.assertBlocked(replace(self.assay, calibrant=calibrant),
                                   "calibrant per-scan area is outside the numeric range", (self.areas[0], None, None))

    def test_summed_rational_area_divides_before_rounding(self):
        area = Fraction(5, 2**1076)
        sample = replace(self.assay.sample, area=area, scans=2)
        calibrant = replace(self.assay.calibrant, area=area, scans=2.0)
        result = calculate_external_calibration(replace(self.assay, sample=sample, calibrant=calibrant))
        self.assertEqual((result.reasons, result.J_a, result.J_s, result.R), ((CALCULATION_ONLY,), 5e-324, 5e-324, 1))


class ExternalCalibrationLinkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assay, _ = load_assay_case()
        cls.original = calculate_external_calibration(cls.assay)
        cls.prerequisite = Prerequisite("pair", "pass", "Documented")
        cls.starting = replace(cls.assay, **{
            role: replace(getattr(cls.assay, role), area=1, scans=None, receiver_gain=None, reader_scale=None,
                          pulse_sequence=None, temperature=None, storage_scale=None, calibrated_flip_angle_deg=None)
            for role in ("sample", "calibrant")
        }, prerequisites=(cls.prerequisite,))
        raw, clean, cls.quiet = clean_regions()
        cls.inputs = {}
        for role in ("sample", "calibrant"):
            signals = getattr(cls.assay, f"{role}_signals")
            run = getattr(cls.assay, role)
            region = replace(clean, nominal_area=run.area / len(signals), run=run.run_id)
            cls.inputs[f"{role}_raw"] = replace(raw, source=Path(run.run_id), scans=run.scans,
                                                receiver_gain=run.receiver_gain, pulse_program=run.pulse_sequence,
                                                temperature_k=298.15, intensity_scale=1.0,
                                                relaxation_delay_s=run.relaxation_delay_s, storage_scale=run.storage_scale)
            cls.inputs[f"{role}_regions"] = tuple((signal, region) for signal in signals)
            cls.inputs[f"{role}_acquisition"] = acquisition_checks(cls.inputs[f"{role}_raw"], signals)

    def link(self, assay=None, **changes):
        return link_external_checks(self.starting if assay is None else assay, **{**self.inputs, **changes})

    def assertLinkBlocked(self, assay, reason, *, keep_areas=True):
        result = calculate_external_calibration(assay)
        self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
        self.assertIn(reason, result.reasons)
        if keep_areas:
            self.assertEqual((result.J_a, result.J_s, result.R), (self.original.J_a, self.original.J_s, self.original.R))

    def test_core_has_no_status_for_entered_linked_or_relinked_areas(self):
        linked = self.link()
        for record in (self.assay, linked, self.link(linked)):
            calculation = _calculate_external_calibration(record)
            self.assertFalse(hasattr(calculation, "status"))
            self.assertEqual((calculation.value, calculation.reasons), (self.original.value, ()))
            self.assertEqual(calculate_external_calibration(record).status, "cannot be assessed")

    def test_linked_runs_fill_headers(self):
        linked = self.link(sample_regions=iter(self.inputs["sample_regions"]))
        recorded = {name: getattr(self.inputs["sample_raw"], name) for name in RECORDED}
        self.assertEqual((linked.sample, linked.calibrant),
                         (replace(self.assay.sample, **recorded), replace(self.assay.calibrant, **recorded)))
        self.assertEqual(linked.prerequisites[0], self.prerequisite)
        self.assertEqual(len(linked.prerequisites), 9)
        self.assertTrue(all(check.status == "pass" for check in linked.prerequisites))
        self.assertEqual(calculate_external_calibration(linked), self.original)

    def test_link_external_checks_handles_missing_run_records(self):
        for starting in (self.starting, self.link()):
            for role, area in (("sample", "J_s"), ("calibrant", "J_a")):
                with self.subTest(linked=starting is not self.starting, role=role):
                    linked = self.link(replace(starting, **{role: None}))
                    missing = f"{role} run: Missing run record"
                    for assay in (linked, self.link(linked)):
                        result = calculate_external_calibration(assay)
                        self.assertEqual(result.reasons.count(missing), 1)
                        self.assertEqual(getattr(result, area), getattr(self.original, area))
                        self.assertIsNone(result.R)
                        self.assertIsNone(purity_report(result)[role]["header"])

    def test_removed_run_preserves_surviving_region_unit_warning(self):
        (signal, region), = self.inputs["calibrant_regions"]
        linked = self.link(calibrant_regions=((signal, replace(region, area=replace(region.area, unit=None))),))
        warning = "linked area units: Region units are missing"
        self.assertIn(f"{warning} (I149 sample region 1 = 'intensity*Hz', {signal} calibrant region 1 = None)",
                      calculate_external_calibration(linked).reasons)
        result = calculate_external_calibration(replace(linked, sample=None))
        self.assertIn(f"{warning} ({signal} calibrant region 1 = None)", result.reasons)
        self.assertEqual(result.J_s, self.original.J_s)
        result = calculate_external_calibration(replace(linked, calibrant=None))
        self.assertFalse(any(reason.startswith(warning) for reason in result.reasons))

    def test_removed_run_respects_current_prerequisites_and_order(self):
        linked = self.link()
        recovery, snr, resolution = (check for check in linked.prerequisites if check.name.startswith("calibrant "))
        sample = next(check for check in linked.prerequisites if check.name.startswith("sample "))
        flag = replace(recovery, status="flag", reason="Reviewed recovery evidence")
        current = (resolution, flag, sample, snr, self.prerequisite)
        result = calculate_external_calibration(replace(linked, sample=None, prerequisites=current))
        self.assertIn(flag.message, result.reasons)
        missing = Prerequisite("sample run", "cannot be assessed", "Missing run record")
        self.assertEqual(purity_report(result)["prerequisites"],
                         [vars(check) for check in (flag, self.prerequisite, resolution, snr, missing)])

    def test_two_linked_files_are_not_transfer_evidence(self):
        self.assertLinkBlocked(self.link(replace(self.starting, pair=None)), "pair record is missing")
        self.assertLinkBlocked(self.link(replace(self.starting, S=None)),
                               "equal nominal receiver gains are not scale evidence")
        unknown = {f"{role}_raw": replace(self.inputs[f"{role}_raw"], intensity_scale=None)
                   for role in ("sample", "calibrant")}
        self.assertLinkBlocked(self.link(**unknown), "sample reader_scale is missing", keep_areas=False)
        fixtures = Path(__file__).parent / "fixtures"
        for raw in (read_jeol(fixtures / "jeol" / "nishizaki_fid02_ma_dmso2_d20s_ns32.jdf"),
                    read_bruker(fixtures / "bruker" / "zenodo17083860_833-1" / "10")):
            self.assertEqual(raw.intensity_scale, 1.0)

    def test_relink_drops_previous_reader_scale_when_current_reader_is_unknown(self):
        linked = self.link()
        raw = replace(self.inputs["sample_raw"], intensity_scale=None)
        self.assertLinkBlocked(self.link(linked, sample_raw=raw), "sample reader_scale is missing", keep_areas=False)
        documented = replace(linked, sample=replace(linked.sample, reader_scale_source="Independent scale log"))
        self.assertEqual(calculate_external_calibration(self.link(documented, sample_raw=raw)), self.original)

    def test_storage_scale_evidence_must_cover_both_runs(self):
        fixture = Path(__file__).parent / "fixtures" / "bruker" / "zenodo17083860_833-1" / "10"
        raw = read_bruker(fixture)
        self.assertEqual(raw.storage_scale, "NC=-6")
        old = (raw.storage_scale, raw.storage_scale)
        assay = replace(self.starting, S_storage_scales=old)
        sample_raw = replace(self.inputs["sample_raw"], storage_scale="NC=-5")
        calibrant_raw = replace(self.inputs["calibrant_raw"], storage_scale=raw.storage_scale)
        linked = self.link(assay, sample_raw=sample_raw, calibrant_raw=calibrant_raw)
        current = ("NC=-5", "NC=-6")
        self.assertLinkBlocked(linked, f"S_storage_scales {old!r} do not cover the run storage scales {current!r}")
        self.assertEqual(calculate_external_calibration(replace(linked, S_storage_scales=current)), self.original)
        unknown = self.link(linked, sample_raw=replace(sample_raw, storage_scale=None), calibrant_raw=calibrant_raw)
        self.assertLinkBlocked(unknown, "sample storage_scale is missing")

    def test_linked_flip_angle_mismatch_needs_correction(self):
        recovery, *other = self.inputs["sample_acquisition"]
        raw = replace(self.inputs["sample_raw"], pulse_width_us=10.0)
        signal = SignalEvidence(recovery.signal, 2.0, "T1 log")
        for calibration in (RecoveryEvidence(calibrated_flip_angle_deg=30.0, calibration_source="Pulse log"),
                            RecoveryEvidence(calibrated_90_us=30.0, calibration_source="Pulse log")):
            with self.subTest(calibration=calibration):
                check = check_recovery(raw, signal, calibration)
                self.assertEqual(check.status, "pass")
                linked = self.link(sample_acquisition=(check, *other))
                self.assertLinkBlocked(linked, "calibrated_flip_angle_deg differs without K evidence")
                corrected = replace(linked, K_covers=("calibrated_flip_angle_deg",))
                self.assertEqual(calculate_external_calibration(corrected), self.original)
                self.assertLinkBlocked(replace(corrected, K=replace(corrected.K, source="NOT FOUND")),
                                       "K.source is missing")

    def test_linked_run_sums_two_disjoint_regions(self):
        raw = self.inputs["sample_raw"]
        names = ("anchor", "singlet")
        bounds = tuple(ppm_bounds((center - 50, center + 50)) for center in (-900, 0))
        spectrum = process(raw, QUANTITATIVE_SETTINGS, regions_ppm=bounds)
        regions = []
        for name, span, center in zip(names, bounds, ppm_bounds((-900, 0))):
            evidence = ComponentEvidence((Line(center, 1.0, 1.0),), "fixture",
                                         satellite_convention="absent", satellite_source="fixture")
            region = check_region(spectrum, span, evidence, raw=raw, noise_bounds_ppm=ppm_bounds(NOISE_BOUNDS_HZ))
            self.assertGreater(region.nominal_area, 0)
            regions.append((name, region))
        expected = sum(region.nominal_area for _, region in regions)
        linked = self.link(replace(self.starting, sample_signals=names), sample_regions=iter(regions),
                           sample_acquisition=acquisition_checks(raw, names))
        self.assertEqual(linked.sample.area, expected)
        result = calculate_external_calibration(linked)
        self.assertEqual((result.J_a, result.reasons), (expected / raw.scans, (CALCULATION_ONLY,)))

    def test_acquisition_ignores_other_runs_unselected_signal(self):
        unrelated = replace(self.inputs["sample_acquisition"][-1], status="flag")
        checks = (*self.inputs["calibrant_acquisition"], unrelated)
        self.assertEqual(self.link(calibrant_acquisition=checks), self.link())

    def test_required_regions_and_checks(self):
        for role in ("sample", "calibrant"):
            key = f"{role}_regions"
            regions = self.inputs[key]
            for index, (signal, _) in enumerate(regions):
                with self.subTest(role=role, signal=signal):
                    linked = self.link(**{key: regions[:index] + regions[index + 1:]})
                    self.assertLinkBlocked(linked, f"{signal} {role} region: Missing region", keep_areas=False)
            key = f"{role}_acquisition"
            checks = self.inputs[key]
            for index, check in enumerate(checks):
                with self.subTest(role=role, signal=check.signal, check=check.check):
                    linked = self.link(**{key: checks[:index] + checks[index + 1:]})
                    self.assertLinkBlocked(linked, f"{role} {check.signal} {check.check}: Missing {check.check} check")

    def test_blank_signal_names_block(self):
        (_, region), = self.inputs["calibrant_regions"]
        checks = tuple(replace(check, signal="") for check in self.inputs["calibrant_acquisition"])
        linked = self.link(calibrant_regions=(("", region),), calibrant_acquisition=checks)
        self.assertLinkBlocked(linked, "calibrant region 1: Region signal name is blank")
        self.assertLinkBlocked(linked, "calibrant recovery: Check signal name is blank")
        self.assertLinkBlocked(self.link(replace(self.starting, calibrant_signals=())),
                               "calibrant signals: Declare nonblank signal names")
        linked = self.link(replace(self.starting, calibrant_signals=("NF",)), calibrant_regions=(("NF", region),),
                           calibrant_acquisition=tuple(replace(check, signal="NF") for check in checks))
        for reason in ("calibrant region 1: Region signal name is blank", "calibrant recovery: Check signal name is blank",
                       "calibrant signals: Declare nonblank signal names"):
            self.assertLinkBlocked(linked, reason)

    def test_run_id_mismatch_leaves_record_unfilled_and_can_be_corrected(self):
        for role in ("sample", "calibrant"):
            with self.subTest(role=role):
                original = getattr(self.starting, role)
                run = replace(original, run_id="other.jdf", area=None)
                blocked = self.link(replace(self.starting, **{role: run}))
                self.assertEqual(getattr(blocked, role), run)
                reason = f"record 'other.jdf' differs from raw source {original.run_id!r}"
                self.assertLinkBlocked(blocked, f"{role} run_id: {reason}", keep_areas=False)
                baseline = self.link()
                other = "sample" if role == "calibrant" else "calibrant"
                mismatch = Prerequisite(f"{role} run_id", "cannot be assessed", reason)
                self.assertIn(mismatch, blocked.prerequisites)
                self.assertEqual([p for p in blocked.prerequisites if p != mismatch],
                                 [p for p in baseline.prerequisites if other in p.name or p == self.prerequisite])
                corrected = self.link(replace(blocked, **{role: replace(run, run_id=original.run_id)}))
                self.assertEqual(corrected, baseline)
                self.assertEqual(calculate_external_calibration(corrected), self.original)

    def test_linked_tiny_rational_sum_blocks_with_numeric_range(self):
        area = Fraction(1, 10**400)
        (signal, region), = self.inputs["sample_regions"]
        regions = ((signal, replace(region, nominal_area=area)),)
        linked = self.link(sample_regions=regions)
        self.assertEqual(linked.sample.area, area)
        result = calculate_external_calibration(linked)
        self.assertEqual((result.status, result.J_a, result.J_s, result.R),
                         ("cannot be assessed", None, self.original.J_s, None))
        self.assertIn("sample per-scan area is outside the numeric range", result.reasons)

    def test_foreign_flag_stays_with_its_run(self):
        (signal, region), = self.inputs["calibrant_regions"]
        foreign = replace(region, status="flag", flags=(("baseline", "distorted"),), run="other.jdf")
        linked = self.link(calibrant_regions=((signal, foreign),))
        name = f"{signal} calibrant region 1"
        mismatch = f"Run mismatch: expected {self.assay.calibrant.run_id!r}, got 'other.jdf'"
        self.assertLinkBlocked(linked, f"{name}: {mismatch}")
        self.assertEqual([p for p in linked.prerequisites if p.name == name],
                         [Prerequisite(name, "cannot be assessed", mismatch)])
        self.assertIn(Prerequisite(f"other.jdf {name}", "flag", "baseline: distorted"), linked.prerequisites)

    def test_linked_blockers_and_restoration(self):
        (signal, region), = self.inputs["calibrant_regions"]
        run = self.assay.calibrant.run_id
        recovery, snr, resolution = self.inputs["calibrant_acquisition"]
        region_name = f"{signal} calibrant region 1"
        cases = (
            ("calibrant_regions", ((signal, replace(region, status="flag", flags=(("baseline", "distorted"),),
                                                    gaps=(("capture", "no evidence"),))),),
             f"{region_name}: baseline: distorted. capture: no evidence"),
            ("calibrant_regions", ((signal, replace(region, run=None)),),
             f"{region_name}: region run unknown"),
            ("calibrant_regions", ((signal, replace(region, area=replace(region.area, unit="intensity*ppm"))),),
             "linked area units: Region units differ (I149 sample region 1 = 'intensity*Hz', "
             "H146 calibrant region 1 = 'intensity*ppm')"),
            ("calibrant_acquisition", (replace(recovery, run="other.jdf"), snr, resolution),
             f"calibrant H146 recovery: Run mismatch: expected {run!r}, got 'other.jdf'"),
            ("calibrant_acquisition", (recovery, replace(snr, status="flag", reason="Low S/N"), resolution),
             "calibrant H146 snr: Low S/N"),
            ("calibrant_acquisition", (replace(recovery, inputs=()), snr, resolution),
             "calibrant H146 recovery: t1_source is missing. calibration_source is missing"),
        )
        for field, value, reason in cases:
            with self.subTest(reason=reason):
                self.assertLinkBlocked(self.link(**{field: value}), reason)
                self.assertEqual(calculate_external_calibration(self.link(**{field: self.inputs[field]})), self.original)

    def test_missing_facts_block_only_dependent_results(self):
        (signal, region), = self.inputs["calibrant_regions"]
        name, check = f"{signal} calibrant region 1", f"calibrant {signal}"
        units = "Region units are missing (I149 sample region 1 = 'intensity*Hz', H146 calibrant region 1 = None)"

        def checks(**changes):
            return {"calibrant_acquisition": acquisition_checks(self.inputs["calibrant_raw"], (signal,), **changes)}

        def regions(base=region, **changes):
            return {"calibrant_regions": ((signal, replace(base, nominal_area=region.nominal_area, **changes)),)}

        for narrow in ((), ("I149",)):
            base = {"sample_acquisition": acquisition_checks(self.inputs["sample_raw"], ("I149",), narrow=narrow),
                    **checks()}
            baseline = self.link(**base)
            expected = calculate_external_calibration(baseline)
            if narrow:
                self.assertIn("sample I149 resolution: Outside the pilot limit", expected.reasons)
            else:
                self.assertEqual(expected, self.original)
            cases = (
                ({"calibrant_acquisition": checks()["calibrant_acquisition"][1:]},
                 [(f"{check} recovery", "Missing recovery check")]),
                ({**regions(self.quiet, run=region.run), **checks(noise=None)},
                 [(name, NOISE_GAP), (f"{check} snr", "noise_sd is missing")]),
                (checks(no_t1_source=(signal,)), [(f"{check} recovery", "t1_source is missing")]),
                (checks(calibration_source=None), [(f"{check} recovery", "calibration_source is missing")]),
                *((regions(area=replace(region.area, unit=unit)), [("linked area units", units)]) for unit in ("", None, "NF")),
                (regions(run=None), [(name, "region run unknown")]),
            )
            for changes, added in cases:
                with self.subTest(narrow=narrow, fact=added[0]):
                    linked = self.link(**{**base, **changes})
                    added = [Prerequisite(name, "cannot be assessed", reason) for name, reason in added]
                    self.assertEqual([p for p in linked.prerequisites if p not in baseline.prerequisites], added)
                    names = {p.name for p in added}
                    self.assertLessEqual({p for p in baseline.prerequisites if p.name not in names},
                                         set(linked.prerequisites))
                    expected_angle = (None if any(p.reason == "Missing recovery check" for p in added)
                                      else baseline.calibrant.calibrated_flip_angle_deg)
                    self.assertEqual(linked.sample, baseline.sample)
                    self.assertEqual(linked.calibrant, replace(baseline.calibrant, calibrated_flip_angle_deg=expected_angle))
                    self.assertLinkBlocked(linked, f"{added[0].name}: {added[0].reason}")
                    self.assertEqual(calculate_external_calibration(self.link(**base)), expected)
            missing_id = self.link(replace(self.starting, calibrant=replace(self.starting.calibrant, run_id=None)), **base)
            self.assertIsNone(missing_id.calibrant.scans)
            self.assertLinkBlocked(missing_id, "calibrant run_id is missing", keep_areas=False)

    def test_header_conflicts_and_restoration(self):
        for role in ("sample", "calibrant"):
            run = getattr(self.assay, role)
            for name, record in (("scans", 64), ("receiver_gain", 31), ("pulse_sequence", "zg30"), ("temperature", 26),
                                 ("reader_scale", 2)):
                with self.subTest(role=role, field=name):
                    evidence = {"reader_scale_source": "Independent scale log"} if name == "reader_scale" else {}
                    changed = replace(self.assay, **{role: replace(run, **{name: record}, **evidence)})
                    linked = self.link(changed)
                    self.assertEqual(getattr(getattr(linked, role), name), record)
                    header = {"temperature": 25.0, "reader_scale": 1.0}.get(name, getattr(run, name))
                    result = calculate_external_calibration(linked)
                    self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
                    self.assertIn(f"{role} {name}: record {record} differs from header {header}", result.reasons)
                    self.assertEqual(calculate_external_calibration(self.link(replace(changed, **{role: run}))),
                                     self.original)

    def test_fractional_header_temperature(self):
        for kelvin, celsius in ((295.34999999999997, 22.2), (298.5974, 25.4474)):
            with self.subTest(kelvin=kelvin):
                assay = replace(self.starting, sample=replace(self.starting.sample, temperature=celsius))
                raws = {name: replace(self.inputs[name], temperature_k=kelvin)
                        for name in ("sample_raw", "calibrant_raw")}
                linked = self.link(assay, **raws)
                self.assertEqual(linked.calibrant.temperature, celsius)
                self.assertEqual(calculate_external_calibration(linked), self.original)
        assay = replace(self.starting, sample=replace(self.starting.sample, temperature=25.0))
        raw = replace(self.inputs["sample_raw"], temperature_k=298.5974)
        self.assertLinkBlocked(self.link(assay, sample_raw=raw),
                               "sample temperature: record 25.0 differs from header 25.4474")

    def test_missing_area_mode(self):
        assay = replace(self.starting, calibrant=replace(self.starting.calibrant, area_mode=None))
        result = calculate_external_calibration(self.link(assay))
        self.assertEqual((result.status, result.value), ("cannot be assessed", None))
        self.assertIn("calibrant area_mode is missing", result.reasons)
