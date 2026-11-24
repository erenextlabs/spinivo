import json
import unittest
from contextlib import chdir
from copy import deepcopy
from dataclasses import asdict, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from spinivo import RecordError, assess_record, load_record, parse_record
from spinivo.assay_record import _json_object
from spinivo.assessment_report import assessment_report, render_json, render_text
from spinivo.external_calibration import PairRecord
from spinivo.internal_standard import CALCULATION_ONLY, calculate_internal_standard
from tests import test_pipeline as spectral
from tests.integral_fixtures import NOISE_BOUNDS_HZ, ppm_bounds
from tests.test_internal_standard import load_assay_case
from tools.compare_references import load_references

RECORDS = Path(__file__).parent / "fixtures" / "records"
JEOL, EXTERNAL, ENTERED = "jeol_internal_standard", "external_calibration", "fid06_calculation_only"
SHIFT_REASON = "Declare shift_ppm with its source"


def record_data(name, path=(), value=None):
    """Return a fixture record, with the entry at path set to value when a path is given."""
    data = target = json.loads((RECORDS / f"{name}.json").read_text(encoding="utf-8"))
    for key in path[:-1]:
        target = target[key]
    if path:
        target[path[-1]] = value
    return data


class AssayRecordTests(unittest.TestCase):
    def test_spectral_and_calculation_only_records_parse(self):
        internal, external, entered = (load_record(RECORDS / f"{name}.json") for name in (JEOL, EXTERNAL, ENTERED))
        self.assertEqual((internal.runs["run"]["run_id"], internal.assay.analyte_tube_id, internal.assay.standard_tube_id),
                         ("FID02-Pattern1-0619", "FID02-Pattern1-tube", "FID02-Pattern1-tube"))
        sample, calibrant = (str(external.inputs[f"{role}_fid"]) for role in ("sample", "calibrant"))
        pair = PairRecord("record", sample, "sample-prep", calibrant, "calibrant-prep",
                          str(RECORDS / "illustrative_pair_record.txt"))
        self.assertEqual((external.assay.sample.run_id, external.assay.pair), (sample, pair))
        self.assertEqual((entered.mode, entered.inputs, entered.assay.analyte_signals),
                         ("calculation_only", {}, ("H5", "I5")))

    def test_layout_errors_name_the_field(self):
        bounds, name = ("signals", "standard", 0, "bounds_ppm"), ("signals", "standard", 0, "name")
        for record, path, value, message in (
            (JEOL, ("mode",), "internal", "mode: must be one of internal_standard, external_calibration, calculation_only"),
            (JEOL, ("run", "fiel"), "x.jdf", "run.fiel: unknown field"),
            (JEOL, ("signals", "analyte", 0, "area"), 1.0, "signals.analyte[0].area: unknown field"),
            (ENTERED, ("run", "file"), "x.jdf", "run.file: unknown field"),
            (JEOL, ("N_a",), "2", "N_a: must be a number or null"),
            (JEOL, bounds, [2.9, 3.0, 3.01], "signals.standard[0].bounds_ppm: must be a list of 2 items or null"),
            (JEOL, name, "MA", "signals.standard[0].name: repeats 'MA' from signals.analyte[0].name"),
            (EXTERNAL, ("pair",), "pairs/fid02.pdf", "pair: must be an object or null"),
            (EXTERNAL, ("scale", "direction"), "calibrant/sample", "scale.direction: must be 'sample/calibrant'"),
            (JEOL, ("signals", "analyte", 0, "evidence"), {"target_lines": [None]},
             "signals.analyte[0].evidence.target_lines[0]: must be an object"),
        ):
            with self.subTest(message=message):
                with self.assertRaises(RecordError) as caught:
                    parse_record(record_data(record, path, value), RECORDS)
                self.assertEqual(str(caught.exception), message)
        path = Path(self.enterContext(TemporaryDirectory())) / "record.json"
        path.write_text("{", encoding="utf-8")
        with self.assertRaisesRegex(RecordError, "^record: invalid JSON: Expecting property name"):
            load_record(path)

    def test_numbers_outside_float_range_name_the_field(self):
        for value in (float("nan"), float("inf"), 10 ** 400):
            with self.subTest(value=value), self.assertRaisesRegex(RecordError, r"^N_a: .*finite"):
                parse_record(record_data(JEOL, ("N_a",), value), RECORDS)
        path = Path(self.enterContext(TemporaryDirectory())) / "record.json"
        data = json.dumps(record_data(JEOL, ("N_a",), "_number_"))
        for token in ("NaN", "Infinity", "1e400", "1e-400"):
            path.write_text(data.replace('"_number_"', token), encoding="utf-8")
            with self.subTest(token=token), self.assertRaisesRegex(RecordError, r"^N_a: .*finite"):
                load_record(path)

    def test_repeated_members_name_the_field(self):
        path = Path(self.enterContext(TemporaryDirectory())) / "record.json"
        data = json.dumps(record_data(JEOL, ("N_a",), 2))
        for old, new, field in (('"N_a": 2', '"N_a": NaN, "N_a": 2', "N_a"),
                                ('"run": {', '"run": {"file": "x.jdf", ', "run.file")):
            path.write_text(data.replace(old, new, 1), encoding="utf-8")
            with self.subTest(field=field), self.assertRaisesRegex(RecordError, rf"^{field}: repeated field$"):
                load_record(path)

    def test_only_the_first_repeated_member_is_retained(self):
        first = ("schema", "spinivo.assay/1"), ("mode", "internal_standard")
        members = _json_object([*first, *reversed(first), *([first[1]] * 500)])
        self.assertEqual(members.repeated, "mode")
        with self.assertRaisesRegex(RecordError, r"^mode: repeated field$"):
            parse_record(members, RECORDS)

    def test_load_binds_paths_even_after_changing_directory(self):
        with chdir(RECORDS.parent):
            record = load_record(Path("records") / f"{EXTERNAL}.json")
        self.assertEqual(record.assay.pair.source, str(RECORDS / "illustrative_pair_record.txt"))
        for role in ("sample", "calibrant"):
            self.assertEqual(record.inputs[f"{role}_fid"], RECORDS / record_data(EXTERNAL)[role]["file"])

    def test_dot_is_the_record_folder(self):
        internal = parse_record(record_data(JEOL, ("run", "file"), "."), RECORDS)
        self.assertEqual(internal.inputs["fid"], RECORDS)

    def test_calculation_only_reproduces_the_fixture_arithmetic(self):
        result = assess_record(load_record(RECORDS / f"{ENTERED}.json"))
        assay, expected = load_assay_case()
        direct = calculate_internal_standard(assay)
        self.assertEqual((result.status, result.reasons), ("cannot be assessed", (CALCULATION_ONLY,)))
        self.assertEqual((result.R, result.value, result.budget), (direct.R, direct.value, direct.budget))
        self.assertAlmostEqual(result.value, expected["purity"], delta=1e-9)

    def test_jeol_record_runs_through_the_pipeline(self):
        result = assess_record(load_record(RECORDS / f"{JEOL}.json")).result
        ratio = 1 / load_references()["fixtures"]["jeol"]["observed_ratio"]
        self.assertEqual((result.status, result.value), ("cannot be assessed", None))
        self.assertIn("MA/DMSO2 capture_pair: Missing or unusable analyte capture. Missing or unusable standard capture",
                      result.reasons)
        self.assertAlmostEqual(result.R, ratio, delta=1e-9 * ratio)
        self.assertFalse(any("recovery" in reason or "reference" in reason for reason in result.reasons), result.reasons)

    def test_missing_facts_load_and_cannot_be_assessed(self):
        data = record_data(JEOL, ("run", "reference", "source"), None)
        data["signals"]["analyte"][0]["t1_s"] = data["signals"]["standard"][0]["t1_source"] = None
        del data["run"]["recovery"]["extra_recovery_s"]
        result = assess_record(parse_record(data, RECORDS)).result
        self.assertEqual((result.status, result.R is None), ("cannot be assessed", False))
        for reason in ("MA recovery: t1_s is missing", "DMSO2 recovery: t1_source is missing",
                       f"reference: {SHIFT_REASON}", "extra_recovery_s: Must be finite and nonnegative"):
            self.assertIn(reason, result.reasons)
        data = record_data(EXTERNAL, ("pair",), None)
        del data["quantities"]["c_s"]
        del data["calibrant"]["reference"]["shift_ppm"]
        missing, complete = (assess_record(parse_record(item, RECORDS)).result for item in (data, record_data(EXTERNAL)))
        self.assertEqual(missing.status, "cannot be assessed")
        for reason in ("pair record is missing", "c_s is missing", f"calibrant reference: {SHIFT_REASON}"):
            self.assertIn(reason, missing.reasons)
            self.assertNotIn(reason, complete.reasons)


class SpectralRecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        internal, external = spectral.PipelineTests, spectral.ExternalPipelineTests
        internal.setUpClass()
        external.setUpClass()
        cls.raws = {"internal.jdf": internal.raw, **{f"{role}.jdf": raw for role, raw in external.raws.items()}}
        cls.internal = record_data(ENTERED)
        cls.internal["mode"] = "internal_standard"
        cls.external = record_data(EXTERNAL)
        for data, method, runs in ((cls.internal, internal, (("run", "internal"),)),
                                   (cls.external, external, (("sample", "sample"), ("calibrant", "calibrant")))):
            for run, stem in runs:
                data[run].update(file=f"{stem}.jdf", noise_bounds_ppm=ppm_bounds(NOISE_BOUNDS_HZ),
                                 reference={"shift_ppm": 0, "source": "Synthetic carrier"},
                                 recovery=asdict(spectral.RECOVERY))
            for role in data["signals"]:
                names = getattr(method.assay, f"{role}_signals")
                signals = method.signals[role] if method is external else (method.signals[name] for name in names)
                data["signals"][role] = [asdict(signal) for signal in signals]
        fields = ("area_mode", "instrument", "probe", "decoupling", "excitation", "solvent", "tube", "fill")
        for role in ("sample", "calibrant"):
            run = getattr(external.assay, role)
            cls.external[role]["conditions"] = {name: getattr(run, name) for name in fields}
        cls.external["scale"].update(gains=[30, 30], storage_scales=["JEOL", "JEOL"])
        cls.expected = {"internal_standard": internal.expected["R"], "external_calibration": external.expected["R"]}

    def assessment(self, data):
        # Normalize test dataclass tuples into the JSON arrays a real record supplies.
        record = parse_record(json.loads(json.dumps(data)), RECORDS)
        with patch("spinivo.pipeline.read", side_effect=lambda path: replace(self.raws[path.name], source=path)):
            return assess_record(record)

    def assess(self, data):
        return self.assessment(data).result

    def test_report_shows_each_runs_record_ids_next_to_its_file(self):
        for original, labels in ((self.internal, ("run",)), (self.external, ("sample", "calibrant"))):
            for missing in (False, True):
                data = deepcopy(original)
                ids = {label: (None, None) if missing else (data[label]["run_id"], data[label]["preparation_id"])
                       for label in labels}
                for label in labels if missing else ():
                    del data[label]["run_id"]
                    data[label]["preparation_id"] = "n/a"
                with self.subTest(mode=data["mode"], missing=missing):
                    assessment = self.assessment(data)
                    report = json.loads(render_json(assessment_report(assessment)))
                    text = render_text(report)
                    self.assertEqual(report["status"] == "pass", not missing, report["reasons"])
                    self.assertEqual(len(set(sum(ids.values(), ()))), 1 if missing else 2 * len(labels))
                    for label, (run_id, preparation_id) in ids.items():
                        run = report["runs"][label]
                        self.assertEqual((run["source"], run["run_id"], run["preparation_id"]),
                                         (str(RECORDS / data[label]["file"]), run_id, preparation_id))
                        self.assertIn(f'  Run {label}: {run["source"]}, run_id = {run_id or "missing"}, '
                                      f'preparation_id = {preparation_id or "missing"}\n', text)
                    if missing and labels != ("run",):
                        self.assertIn("sample run_id: is missing", report["reasons"])

    def test_entered_run_id_is_the_records_own_id(self):
        data = deepcopy(self.external)
        data["sample"] = {"file": data["sample"]["file"]}
        report = assessment_report(self.assessment(data))["method_report"]
        sample, calibrant = report["sample"], report["calibrant"]
        self.assertEqual((sample["run_id"], sample["entered"]), (str(RECORDS / data["sample"]["file"]), None))
        self.assertEqual((calibrant["run_id"], calibrant["entered"]["run_id"]),
                         (str(RECORDS / data["calibrant"]["file"]), "calibrant-run"))
        self.assertIn("sample run_id: is missing", report["reasons"])
        del data["calibrant"]["run_id"]
        entered = assessment_report(self.assessment(data))["method_report"]["calibrant"]["entered"]
        self.assertEqual((entered["preparation_id"], entered.get("run_id")), ("calibrant-prep", None))

    def test_both_modes_need_evidence_and_pass_when_it_is_restored(self):
        for original, section, field, reason in (
            (self.internal, "run", "noise_bounds_ppm", "noise bounds: Declare two finite ppm bounds"),
            (self.external, "scale", "direction", "scale.direction: is missing"),
        ):
            data = deepcopy(original)
            with self.subTest(mode=data["mode"]):
                complete = self.assess(data)
                self.assertEqual(complete.status, "pass", complete.reasons)
                evidence = data[section].pop(field)
                if section == "run":
                    data[section][field] = [None, evidence[1]]
                result = self.assess(data)
                self.assertEqual((result.status, result.value), ("cannot be assessed", None))
                self.assertIn(reason, result.reasons)
                data[section][field] = evidence
                self.assertEqual(self.assess(data), complete)

    def test_pair_run_labels_are_checked_before_using_paths(self):
        data = deepcopy(self.external)
        self.assertEqual(self.assess(data).status, "pass")
        label = data["sample"]["run_id"]
        data["sample"]["run_id"] = data["pair"]["sample_run_id"] = data["calibrant"]["run_id"]
        result = self.assess(data)
        self.assertEqual(result.status, "cannot be assessed")
        self.assertIn("run_id: sample and calibrant run_id must differ", result.reasons)
        del data["sample"]["run_id"]
        data["pair"]["sample_run_id"] = str(RECORDS / data["sample"]["file"])
        result = self.assess(data)
        self.assertEqual(result.status, "cannot be assessed")
        self.assertIn("sample run_id: is missing", result.reasons)
        data["sample"]["run_id"] = label
        self.assertIn("pair.sample_run_id: must match sample.run_id", self.assess(data).reasons)
        data["sample"]["run_id"] = data["pair"]["sample_run_id"] = label
        self.assertEqual(self.assess(data).status, "pass")

    def test_reference_shift_uses_each_runs_declared_ppm_coordinates(self):
        for original, runs in ((self.internal, (("run", ("analyte", "standard"), 0.5),)),
                               (self.external, (("sample", ("sample",), 0.5), ("calibrant", ("calibrant",), -0.5)))):
            data = deepcopy(original)
            for run, roles, shift in runs:
                data[run]["reference"]["shift_ppm"] = shift
                data[run]["noise_bounds_ppm"] = [value + shift for value in data[run]["noise_bounds_ppm"]]
                for role in roles:
                    for signal in data["signals"][role]:
                        signal["bounds_ppm"] = [value + shift for value in signal["bounds_ppm"]]
                        for line in signal["evidence"]["target_lines"]:
                            line["ppm"] += shift
            with self.subTest(mode=data["mode"]):
                result = self.assess(data)
                self.assertEqual(result.status, "pass", result.reasons)
                expected = self.expected[data["mode"]]
                self.assertAlmostEqual(result.R, expected, delta=1e-4 * expected)

    def test_null_component_linewidth_loads_and_blocks(self):
        data = deepcopy(self.internal)
        data["signals"]["analyte"][0]["evidence"]["target_lines"][0]["fwhm_hz"] = None
        result = self.assess(data)
        self.assertEqual((result.status, result.value), ("cannot be assessed", None))
        self.assertIn("signals.analyte[0].evidence.target_lines[0].fwhm_hz: is missing", result.reasons)

    def test_pair_source_trailing_slash_is_rejected(self):
        external = spectral.ExternalPipelineTests
        pair = external.assay.pair
        assay = replace(external.assay, pair=replace(pair, source=pair.source + "/"))
        data = deepcopy(self.external)
        data["pair"]["source"] += "/"
        for route, run in (("pipeline", lambda: external.assess(assay=assay).result),
                           ("record", lambda: self.assess(data))):
            with self.subTest(route=route):
                result = run()
                self.assertEqual((result.status, result.value), ("cannot be assessed", None))
                self.assertTrue(any(reason.startswith("pair source file could not be read:")
                                    for reason in result.reasons), result.reasons)

    def test_empty_pair_source_blocks_pipeline_and_record(self):
        external = spectral.ExternalPipelineTests
        source = Path(self.enterContext(TemporaryDirectory())) / "pair.txt"
        data = deepcopy(self.external)
        data["pair"]["source"] = str(source)
        assay = replace(external.assay, pair=replace(external.assay.pair, source=str(source)))
        for contents in (b"", b"Pair evidence"):
            source.write_bytes(contents)
            for route, assess in (("pipeline", lambda: external.assess(assay=assay).result),
                                  ("record", lambda: self.assess(data))):
                with self.subTest(route=route, contents=contents):
                    result = assess()
                    if contents:
                        self.assertEqual((result.status, result.reasons), ("pass", ()))
                    else:
                        self.assertEqual((result.status, result.value, result.budget),
                                         ("cannot be assessed", None, None))
                        self.assertEqual(result.reasons, ("pair source file is empty",))
                        self.assertTrue(all(value is not None for value in (result.J_a, result.J_s, result.c_s)))

    def test_satellites_require_an_explicit_fraction(self):
        for convention, coupling in (("include", 20), ("exclude", 260)):
            with self.subTest(convention=convention):
                data = deepcopy(self.internal)
                evidence = data["signals"]["analyte"][0]["evidence"]
                evidence.update(satellite_convention=convention, coupling_hz=coupling)
                complete = self.assess(data)
                self.assertEqual(complete.status, "pass", complete.reasons)
                fraction = evidence.pop("satellite_fraction")
                for value in ({}, {"satellite_fraction": None}):
                    evidence.update(value)
                    result = self.assess(data)
                    self.assertEqual((result.status, result.value), ("cannot be assessed", None))
                    self.assertIn("signals.analyte[0].evidence.satellite_fraction: is missing", result.reasons)
                    self.assertIn("H5 analyte region 1: satellites: Missing satellite fraction", result.reasons)
                    self.assertFalse(any("satellite convention" in reason for reason in result.reasons), result.reasons)
                evidence["satellite_fraction"] = fraction
                self.assertEqual(self.assess(data), complete)
