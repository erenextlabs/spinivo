import json
import unittest
from contextlib import redirect_stdout
from copy import deepcopy
from dataclasses import replace
from functools import reduce
from io import StringIO
from operator import getitem
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from spinivo import RecordError, assess_record, parse_record
from spinivo.assessment_report import assessment_report
from spinivo.cli import main
from spinivo.integral_components import check_region
from spinivo.measurement import measure_region
from spinivo.processing import process
from tests import test_assay_record as records
from tests.test_cli import write_synthetic_bruker

FWHM, UNIT, RUN = "measured fwhm_hz", "measured area unit", "whole run"
NO_FWHM = "Missing a half-height crossing inside the signal bounds"
SERIES = "series:fixture-K-D"
# Each run: record key, check name prefix and its signals as (role, index, name).
RUNS = {"internal": (("run", "", (("analyte", 0, "H5"), ("analyte", 1, "I5"), ("standard", 0, "G5"))),),
        "external": (("sample", "sample ", (("sample", 0, "I149"),)),
                     ("calibrant", "calibrant ", (("calibrant", 0, "H146"),)))}
AREA = {"sample": "J_a", "calibrant": "J_s"}
CONDITIONS = ("instrument", "probe", "decoupling", "excitation", "solvent", "tube", "fill")
QUANTITIES = {"internal": ("M_a", "M_s", "m_a", "m_s", "w_s"), "certified": ("c_s",),
              "prepared": ("m_s", "w_s", "M_s", "V_s")}


def spectral_cases(method):
    """Facts each run and signal needs: (fact, removal, expected reasons, names of the evidence that may change)."""
    internal, cases = method == "internal", []
    all_regions = [f"{name} {role} region {index + 1}" for _, _, signals in RUNS[method]
                   for role, index, name in signals]
    for key, prefix, signals in RUNS[method]:
        recovery = tuple(f"{prefix}{name} recovery" for _, _, name in signals)
        regions = tuple(f"{name} {owner} region" for owner, _, name in signals)
        areas = ("analyte_areas", "standard_areas") if internal else ()
        if internal:
            pair_checks = tuple(f"{first}/{second} capture_pair" for role, _, first in signals if role == "analyte"
                                for owner, _, second in signals if owner == "standard")
            run_dependents = (key, "I_a", "I_s", "R", *areas, *(f"source {name}" for name in areas), *regions,
                              *pair_checks, *(f"{name} {check}" for _, _, name in signals
                                              for check in ("recovery", "snr", "resolution")))
        else:
            run_dependents = (key, *regions, "I149/H146 capture_pair", AREA[key], "R")
        cases += [(f"{owner} T1", ("signals", owner, position, field), f"{prefix}{name} recovery: {field} is missing",
                   (f"{prefix}{name} recovery",))
                  for owner, position, name in signals for field in ("t1_s", "t1_source")]
        for role, index, name in signals:
            region = f"{name} {role} region"
            missing = f"{region} {index + 1}"
            units = ", ".join(f"{entry} = {None if entry == missing else 'intensity*Hz'!r}" for entry in all_regions)
            cases += [
                ("resolution input", (FWHM, role, index), f"{prefix}{name} resolution: fwhm_hz is missing: {NO_FWHM}",
                 (f"{prefix}{name} resolution",)),
                ("area units", (UNIT, role, index), f"linked area units: Region units are missing ({units})",
                 ("linked area units", f"{region} values", *areas)),
            ]
        cases += [
            ("pulse calibration source", (key, "recovery", "calibration_source"),
             f"{recovery[0]}: calibration_source is missing", recovery),
            ("pulse calibration", (key, "recovery", "calibrated_flip_angle_deg"),
             f"{recovery[0]}: calibrated_flip_angle_deg or calibrated_90_us is required", recovery),
            ("run ID", (key, "run_id"), "run_id is missing" if internal else f"{key} run_id: is missing", ()),
            ("preparation ID", (key, "preparation_id"),
             "analyte_tube_id is missing" if internal else f"{key} preparation_id is missing", ()),
            ("whole run", (RUN, key), ("reader: Missing raw data file",) if internal else
             (f"{key} run: Missing run record", f"{key} reader: Missing raw data file"),
             run_dependents),
        ]
    return cases


def assay_cases(variant):
    """Facts of the measurement model: (fact, removal, expected reasons, names of the evidence that may change)."""
    internal = variant == "internal"
    concentration = () if internal else QUANTITIES[variant]
    names = QUANTITIES[variant] + (() if internal else ("M_a", "V_a", "m_a", "K", "D", "S"))
    cases = [
        ("repeatability", ("R_relative_u",), "R_relative_u is missing", ()),
        ("repeatability source", ("R_u_source",), "R_u_source is missing", ("source R", "source R_relative_u")),
        ("coverage factor source", ("k", "source"), "k.source is missing", ("source k",)),
        *(("material u", ("scale", "S", "standard_u") if name == "S" else ("quantities", name, "standard_u"),
           f"{name}.standard_u is missing", ()) for name in names),
        *((f"quantity {name}", ("quantities", name), f"{name} is missing",
           (f"source {name}", "c_s") if name in concentration else (f"source {name}",))
          for name in (("m_a", "m_s") if internal else ("m_a", *concentration))),
    ]
    if internal:
        return cases + [
            ("covariance", ("correlations", 0),
             "missing covariance for m_s, m_a (shared source balance:fixture-FID06)", ()),
            ("covariance source", ("correlations", 0, "source"), "correlation m_s, m_a source is missing", ()),
        ]
    cases += [
        ("pair record", ("pair",), "pair record is missing", ()),
        ("pair source", ("pair", "source"), "pair source is missing", ()),
        ("pair kind", ("pair", "kind"), "pair kind must be record or attestation", ()),
        *(("pair role ID", ("pair", f"{role}_{name}"), f"pair {role}_{name} is missing", ())
          for role in AREA for name in ("run_id", "preparation_id")),
        ("c_s route", ("c_s_route",), "c_s_route must be declared as certified or prepared", ("c_s",)),
        ("covariance", ("covariances", 0), f"missing covariance for K, D (shared source {SERIES})", ()),
        ("covariance source", ("covariances", 0, "source"), "covariance K, D source is missing", ()),
        ("scale S", ("scale", "S"), "S is missing", ("source S",)),
        ("scale S source", ("scale", "S", "source"), "S.source is missing", ("source S",)),
        ("scale gains", ("scale", "gains"), "S_gains is missing", ()),
        ("scale storage scales", ("scale", "storage_scales"), "S_storage_scales is missing", ()),
        ("scale direction", ("scale", "direction"), "scale.direction: is missing", ()),
        *((name, ("quantities", name), f"{name} is missing", (f"source {name}",)) for name in ("K", "D")),
        *((f"{name} source", ("quantities", name, "source"), f"{name}.source is missing", (f"source {name}",))
          for name in ("K", "D")),
        ("transfer K evidence", ("K_covers",), "solvent differs without K evidence", ()),
        *(("transfer", (role, "conditions", name), f"{role} {name} is missing",
           (AREA[role], "R", f"{role} area") if name == "area_mode" else ())
          for role in AREA for name in ("area_mode", *CONDITIONS)),
    ]
    if variant == "prepared":
        cases.append(("covariance", ("correlations", 0),
                      "missing covariance for m_s, m_a (shared source balance:fixture-FID02)", ()))
    return cases


class MissingEvidenceTests(unittest.TestCase):
    """Remove one fact at a time from a passing record: purity is withheld for a named reason, then restored."""

    @classmethod
    def setUpClass(cls):
        spectral = records.SpectralRecordTests
        spectral.setUpClass()
        cls.raws, certified = spectral.raws, deepcopy(spectral.external)
        quantities = certified["quantities"]
        quantities["K"]["source_ids"] = quantities["D"]["source_ids"] = [SERIES]
        certified["calibrant"]["conditions"]["solvent"], certified["K_covers"] = "D2O", ["solvent"]
        case = json.loads((records.RECORDS.parent / "external_calibration_cases.json").read_text(encoding="utf-8"))
        cls.prepared_c_s = case["expected"]["c_s"]
        case = case["assay"]
        prepared = deepcopy(certified)
        quantities = prepared["quantities"]
        del quantities["c_s"]
        quantities.update({name: case[name] for name in QUANTITIES["prepared"]})
        for name in ("m_a", "V_a"):
            quantities[name]["source_ids"] = case[name]["source_ids"]
        prepared.update(c_s_route="prepared", correlations=case["correlations"])
        cls.data, cls.spectra = {"internal": spectral.internal, "certified": certified, "prepared": prepared}, {}
        cls.complete = {variant: cls.assessment(cls.record(data)) for variant, data in cls.data.items()}
        cls.evidence_of = {variant: cls.evidence(assessment) for variant, assessment in cls.complete.items()}

    @classmethod
    def record(cls, data):
        return parse_record(json.loads(json.dumps(data)), records.RECORDS)

    @classmethod
    def assessment(cls, record):
        """Assess with each spectrum processed once per file, settings and regions; the removed facts act later."""
        def processed(raw, settings, **options):
            key = (raw.source, settings, *options.items())
            if key not in cls.spectra:
                cls.spectra[key] = process(raw, settings, **options)
            return cls.spectra[key]

        with (patch("spinivo.pipeline.read", side_effect=lambda path: replace(cls.raws[path.name], source=path)),
              patch("spinivo.pipeline.process", side_effect=processed)):
            return assess_record(record)

    @staticmethod
    def evidence(assessment):
        """Return the report and each check, measured result and source in it, by name."""
        report = assessment_report(assessment)
        method = report["method_report"]
        items = [(check["name"], check) for check in method["prerequisites"]]
        items += [(f"{'/'.join(pair['signals'])} capture_pair values", pair) for pair in report["capture_pairs"]]
        items += [(f"source {name}", (quantity["source"], quantity.get("source_ids")))
                  for name, quantity in method["inputs"].items()]
        items.append(("source k", (method["k"] or {}).get("source")))
        for label, run in report["runs"].items():
            prefix = "" if label == "run" else f"{label} "
            items.append((f"{label} source", run["source"]))
            items += [(f"{prefix}{check['signal']} {check['check']} values", check)
                      for check in run["acquisition"] or ()]
            items += [(f"{region['signal']} {region['role']} region values", region)
                      for region in run["regions"] or ()]
        if "I_a" in method:
            items += [(name, method[name]) for name in ("I_a", "I_s", "R")]
            items += [(name, method["inputs"][name]) for name in ("analyte_areas", "standard_areas")]
        else:
            items += [(name, method[name]) for name in ("J_a", "J_s", "R")]
            items.append(("c_s", (method["c_s"], method["concentration_unit"])))
            items += [(f"{role} area", (method[role]["area"], method[role]["per_scan_area"])) for role in AREA]
        return report, items

    @staticmethod
    def measured_without(data, kind, role, index):
        """Drop one measured input of a declared signal: its FWHM or its region's area unit."""
        signal = data["signals"][role][index]
        if kind == FWHM:
            def measure(*args, **kwargs):
                measured = measure_region(*args, **kwargs)
                return replace(measured, fwhm_hz=None, fwhm_reason=NO_FWHM) if args[1] == signal["name"] else measured
            return patch("spinivo.pipeline.measure_region", side_effect=measure)
        file = (data.get(role) or data["run"])["file"]

        def check(spectrum, bounds, *args, **kwargs):
            region = check_region(spectrum, bounds, *args, **kwargs)
            if bounds == tuple(signal["bounds_ppm"]) and kwargs["raw"].source.endswith(file):
                region = replace(region, area=replace(region.area, unit=None))
            return region
        return patch("spinivo.pipeline.check_region", side_effect=check)

    def assertRemoval(self, variant, path, reasons, dependent, *, exact=False):
        data, complete = deepcopy(self.data[variant]), self.complete[variant]
        if path[0] == RUN:
            record, role = self.record(data), path[1]
            inputs = {**record.inputs, "fid" if role == "run" else f"{role}_fid": None}
            assay = record.assay if role == "run" else replace(record.assay, **{role: None})
            blocked = self.assessment(replace(record, assay=assay, inputs=inputs))
        elif path[0] in (FWHM, UNIT):
            with self.measured_without(data, *path):
                blocked = self.assessment(self.record(data))
        else:
            target = reduce(getitem, path[:-1], data)
            value = target.pop(path[-1])
            blocked = self.assessment(self.record(data))
            if isinstance(target, list):
                target.insert(path[-1], value)
            else:
                target[path[-1]] = value
        result, (report, items) = blocked.result, self.evidence(blocked)
        self.assertEqual((result.status, result.value, result.budget), ("cannot be assessed", None, None))
        self.assertEqual((report["status"], report["purity"]["value"], report["method_report"]["budget"]),
                         ("cannot be assessed", None, []))
        if exact:
            self.assertCountEqual(result.reasons, reasons)
        else:
            for reason in (reasons,) if isinstance(reasons, str) else reasons:
                self.assertIn(reason, result.reasons)
        complete_report, before = self.evidence_of[variant]
        for name, item in before:
            if not any(name == entry or name.startswith(f"{entry} ") for entry in dependent):
                self.assertIn((name, item), items, name)
        restored = self.assessment(self.record(data))
        self.assertEqual(restored, complete)
        self.assertEqual(assessment_report(restored), complete_report)

    def assertRemovals(self, variant, cases):
        for fact, path, reasons, dependent in cases:
            with self.subTest(variant=variant, fact=fact, removed=path):
                self.assertRemoval(variant, path, reasons, dependent)

    def test_complete_records_pass_on_both_concentration_routes(self):
        for variant, assessment in self.complete.items():
            with self.subTest(variant=variant):
                self.assertEqual((assessment.result.status, assessment.result.reasons), ("pass", ()))
        self.assertEqual(self.complete["certified"].result.c_s, 2.7)
        self.assertAlmostEqual(self.complete["prepared"].result.c_s, self.prepared_c_s, delta=1e-9)

    def test_internal_standard_removals(self):
        self.assertRemovals("internal", (*spectral_cases("internal"), *assay_cases("internal")))

    def test_missing_noise_only_blocks_noise_checks(self):
        for variant in self.data:
            for key, prefix, signals in RUNS["internal" if variant == "internal" else "external"]:
                with self.subTest(variant=variant, run=key):
                    regions = [(f"{name} {role} region", index + 1) for role, index, name in signals]
                    snr = [f"{prefix}{name} snr" for _, _, name in signals]
                    reasons = (
                        *(f"{region} {number}: snr: Missing noise window. neighbour: Missing noise window"
                          for region, number in regions),
                        *(f"{name}: noise_sd is missing" for name in snr),
                    )
                    self.assertRemoval(variant, (key, "noise_bounds_ppm"), reasons,
                                       (*snr, *(region for region, _ in regions)), exact=True)

    def test_certified_external_removals(self):
        self.assertRemovals("certified", (*spectral_cases("external"), *assay_cases("certified")))

    def test_prepared_external_removals(self):
        self.assertRemovals("prepared", (*spectral_cases("external"), *assay_cases("prepared")))

    def test_coverage_factor_is_optional(self):
        for variant, complete in self.complete.items():
            with self.subTest(variant=variant):
                data = deepcopy(self.data[variant])
                k = data.pop("k")
                result, budget = self.assessment(self.record(data)).result, complete.result.budget
                self.assertEqual((result.status, result.reasons, result.value), ("pass", (), complete.result.value))
                self.assertEqual((result.budget.standard_u, result.budget.expanded_u, result.budget.expanded_u_reason),
                                 (budget.standard_u, None, "coverage factor k was not supplied"))
                data["k"] = k
                self.assertEqual(self.assessment(self.record(data)), complete)

    def test_record_without_a_run_names_the_missing_file(self):
        for variant, key in (("internal", "run"), ("certified", "sample"), ("prepared", "calibrant")):
            with self.subTest(variant=variant):
                data = deepcopy(self.data[variant])
                del data[key]
                with self.assertRaisesRegex(RecordError, rf"^{key}\.file: must name the data file$"):
                    self.record(data)

    def test_command_line_reports_the_same_gate(self):
        root = Path(self.enterContext(TemporaryDirectory()))
        for name, raw in self.raws.items():
            write_synthetic_bruker(root / name, raw)

        def assess(data):
            path = root / "record.json"
            path.write_text(json.dumps(data), encoding="utf-8")
            with redirect_stdout(StringIO()) as stdout:
                self.assertEqual(main(["assess", str(path), "--format", "json"]), 0)
            return json.loads(stdout.getvalue())

        for variant, path, reason in (("internal", ("run", "run_id"), "run_id is missing"),
                                      ("internal", ("signals", "standard", 0, "t1_source"),
                                       "G5 recovery: t1_source is missing"),
                                      ("certified", ("quantities", "c_s"), "c_s is missing"),
                                      ("prepared", ("pair",), "pair record is missing")):
            with self.subTest(variant=variant, removed=path):
                data = deepcopy(self.data[variant])
                if "pair" in data:
                    data["pair"]["source"] = str(records.RECORDS / data["pair"]["source"])
                    data["scale"]["storage_scales"] = ["NC=0", "NC=0"]
                complete = assess(data)
                self.assertEqual(complete["status"], "pass", complete["reasons"])
                target = reduce(getitem, path[:-1], data)
                value = target.pop(path[-1])
                blocked = assess(data)
                self.assertEqual((blocked["status"], blocked["purity"]["value"]), ("cannot be assessed", None))
                self.assertIn(reason, blocked["reasons"])
                target[path[-1]] = value
                self.assertEqual(assess(data), complete)
