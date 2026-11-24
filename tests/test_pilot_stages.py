import unittest
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from dataclasses import fields, replace
from io import StringIO
from unittest.mock import patch

from spinivo import read
from spinivo.processing import ProcessingSettings, integrate, process
from tests.test_readers import FIXTURES
from tools import pilot_stages
from tools.compare_references import _pilot_paths, load_references, measure_regions


class PilotStageTests(unittest.TestCase):
    def test_each_trial_changes_one_named_setting(self):
        processing = ProcessingSettings(**load_references()["processing"])
        (_, _, baseline), *variants = pilot_stages.stages(processing)
        expected = {
            "scan count": ("scan_normalisation", True),
            "receiver gain": ("gain_normalisation", True),
            "reader scaling": ("undo_reader_scale", True),
            "filter": ("remove_filter_delay", False),
            "axis sampling": ("zero_fill_factor", 1),
            "phase": ("phase_lb_hz", 0.0),
            "region placement": ("region_placement", None),
            "pulse direction": ("pulse_direction", None),
        }
        self.assertEqual({stage for stage, _, _ in variants}, expected.keys())
        self.assertEqual(len({settings for _, _, settings in variants}), 11)
        for stage, detail, settings in variants:
            with self.subTest(detail=detail):
                changed = [f.name for f in fields(settings) if getattr(settings, f.name) != getattr(baseline, f.name)]
                field, value = expected[stage]
                self.assertEqual(changed, [field])
                if value is not None:
                    self.assertEqual(getattr(settings, field), value)

    def test_fixture_report_and_independence_from_purity(self):
        references = load_references()
        pilot = references["pilots"][0]
        fixture = references["fixtures"]["jeol"]
        path = FIXTURES / fixture["path"]
        # Exercise the pilot wiring with the fixture's DMSO2 in place of caffeine.
        pilot["runs"][0]["regions"]["I149"]["bounds_ppm"] = fixture["regions_ppm"]["DMSO2"]
        paths = _pilot_paths(pilot)
        raw = read(path)
        sample = replace(raw, scans=16, intensity_scale=2.0, pulse_width_us=8.0, fid=2 * raw.fid)
        ec = replace(raw, pulse_width_us=4.0)
        measured = measure_regions(path, fixture["regions_ppm"], references)
        area = measured["DMSO2"].area
        control = measured["MA"].area / area
        self.assertAlmostEqual(control, 1 / fixture["observed_ratio"])

        def report(record, catalog):
            with patch.object(pilot_stages, "read", side_effect=(sample, ec)), redirect_stdout(StringIO()) as output:
                pilot_stages.report_pilot(FIXTURES, record, catalog, paths)
            return output.getvalue()

        with patch.object(pilot_stages, "process", wraps=pilot_stages.process) as process_run:
            output = report(pilot, references)
        options = ProcessingSettings(**references["processing"])
        for run in (sample, ec):
            calls = [call for call in process_run.call_args_list if call.args[0].fid is run.fid]
            self.assertEqual(len(calls), 4)
            base, bypass, axis, phase = calls
            for call in (bypass, axis, phase):
                self.assertEqual(call.kwargs, base.kwargs)
            self.assertEqual(base.args[1], options)
            self.assertEqual(bypass.args[0].digital_filter, {"digital_filter": False})
            self.assertEqual(bypass.args[1], options)
            self.assertIs(axis.args[0], run)
            self.assertEqual(axis.args[1], replace(options, zero_fill_factor=1))
            self.assertIs(phase.args[0], run)
            self.assertEqual(phase.args[1], replace(options, phase_lb_hz=0.0))
        altered = deepcopy(pilot)
        altered["purity_reference"] = "invalid"
        altered["runs"][0]["regions"]["I149"]["purity"] = "invalid"
        self.assertEqual(output, report(altered, dict(references, certificate_purity="invalid")))
        self.assertIn(f"sample G149/I149 ratio={control:.8g}", output)
        self.assertIn(f"EC G146/H146 ratio={control:.8g}", output)
        self.assertIn(f"Cross-tube I149/H146 workbook_area_ratio={0.55648986 / 0.35070252:.8g}", output)
        for top, bottom in pilot["ratios"]:
            self.assertRegex(output, rf"{top}/{bottom} Baseline shape \(baseline-shape/1, diagnostic\): "
                                     r"computed, [\d.e+-]+ % relative standard uncertainty")
        self.assertIn(f"baseline | as in compare_references | {2 * area:.8g} | {area:.8g} | 2 | +0.0000%", output)
        self.assertIn(f"divide by NS if stored values are sums | {2 * area / 16:.8g} | {area / 32:.8g} | 4 | +100.0000%", output)
        self.assertIn(f"divide by documented reader multiplier | {area:.8g} | {area:.8g} | 1 | -50.0000%", output)
        self.assertIn(f"nominal width sample/EC, arithmetic only | {2 * area:.8g} | {area:.8g} | 4 | +100.0000%", output)
        self.assertIn(f"nominal width EC/sample, arithmetic only | {2 * area:.8g} | {area:.8g} | 1 | -50.0000%", output)
        self.assertIn("receiver response correction unknown | unknown | unknown | unknown | unknown", output)
        self.assertIn("workbook bounds unknown | unknown | unknown | unknown | unknown", output)
        self.assertIn("filter | bypass JEOL group delay correction (I149: Lineshape autophase", output)
        self.assertIn("axis sampling | zero fill factor 1, same header ppm reference", output)
        self.assertIn("phase | fit on LB0 copy", output)
        windows = tuple(map(tuple, fixture["regions_ppm"].values()))
        trace = process(raw, options, regions_ppm=windows)
        # Five percent of the 0.11 ppm window exceeds one acquired point here.
        expanded = integrate(trace, (2.90, 3.0155)).area
        for role, sample_area, ec_area in (("sample", 2 * expanded, area), ("EC", 2 * area, expanded)):
            row = next(line for line in output.splitlines() if f"move {role} upper edge outward" in line)
            self.assertEqual(row.split(" | ")[2:5],
                             [f"{sample_area:.8g}", f"{ec_area:.8g}", f"{sample_area / ec_area:.8g}"])
        self.assertTrue(output.endswith("  unresolved; EC validation gate stays closed\n"))

    def test_missing_bounds_keep_other_measurements_and_later_pair(self):
        references = load_references()
        references["pilots"] = references["pilots"][:2]
        first, later = references["pilots"]
        fixture = references["fixtures"]["jeol"]
        raw = read(FIXTURES / fixture["path"])
        later["runs"][0]["regions"]["G7"]["bounds_ppm"] = fixture["regions_ppm"]["DMSO2"]
        region = first["runs"][0]["regions"]["I149"]
        verdict = "  unresolved; EC validation gate stays closed\n"
        for fields in ({}, {"bounds_ppm": None}, {"bounds_ppm": [2.90, None]}):
            with self.subTest(fields=fields):
                region.pop("bounds_ppm", None)
                region.update(fields)
                with (patch.object(pilot_stages, "load_references", return_value=references),
                      patch.object(pilot_stages, "read", return_value=raw),
                      redirect_stdout(StringIO()) as output, redirect_stderr(StringIO()) as errors):
                    self.assertEqual(pilot_stages.main([str(FIXTURES)]), 0)
                self.assertEqual(errors.getvalue(), "")
                early, late = output.getvalue().split("FID05 (solvent-transfer challenge)\n")
                self.assertIn("sample G149/I149 ratio=unknown workbook_area_ratio=", early)
                self.assertRegex(early, r"EC G146/H146 ratio=\d+\.\d+")
                self.assertIn("baseline_vs_workbook=unknown", early)
                self.assertIn("pilot bounds (ppm) sample=unknown EC=[2.9, 3.01]", early)
                self.assertIn("I149/H146 Baseline shape (baseline-shape/1, diagnostic): cannot be assessed, "
                              "numerator I149: candidate areas unavailable", early)
                rows = [line.split(" | ") for line in early.splitlines() if line.count(" | ") == 5][1:]
                self.assertEqual(len(rows), 12)
                self.assertNotEqual(rows[0][3], "unknown")
                for _, _, sample, _, ratio, change in rows:
                    self.assertEqual((sample, ratio, change), ("unknown", "unknown", "unknown"))
                self.assertTrue(early.endswith(verdict))
                self.assertIn("Cross-tube G7/G6", late)
                self.assertRegex(late, r"baseline \| as in compare_references \| (\d+\.\d+) \| \1 \| 1 \| \+0\.0000%")
                self.assertTrue(late.endswith(verdict))

    def test_nonpositive_areas_warn_and_keep_signed_values(self):
        references = load_references()
        pilot = references["pilots"][1]
        raw = read(FIXTURES / references["fixtures"]["jeol"]["path"])
        for area, ratio, change in ((-2, "1", "+0.0000%"), (0, "unknown", "unknown")):
            with self.subTest(area=area):
                with (patch.object(pilot_stages, "read", return_value=raw),
                      patch.object(pilot_stages, "integrate") as integrate_region,
                      redirect_stdout(StringIO()) as output):
                    integrate_region.return_value.area = area
                    pilot_stages.report_pilot(FIXTURES, pilot, references, _pilot_paths(pilot))
                self.assertIn("baseline | as in compare_references "
                              "(G7: nonpositive signal area. G6: nonpositive signal area) | "
                              f"{area} | {area} | {ratio} | {change}", output.getvalue())

    def test_preflight_accepts_pinned_roles_including_fid06(self):
        with patch.object(pilot_stages, "report_pilot") as report:
            self.assertEqual(pilot_stages.main([str(FIXTURES)]), 0)
        self.assertEqual([call.args[1]["name"] for call in report.call_args_list], ["FID02", "FID05", "FID06"])
        self.assertIn("BTMSBandINH", report.call_args.args[3][0])

    def test_preflight_refuses_unpinned_pair_before_reading(self):
        references = load_references()
        references["pilots"][1]["runs"][0]["path"] = references["fixtures"]["jeol"]["path"]
        with (patch.object(pilot_stages, "load_references", return_value=references),
              patch.object(pilot_stages, "read") as read_run,
              redirect_stdout(StringIO()) as output, redirect_stderr(StringIO()) as errors,
              self.assertRaises(SystemExit) as stopped):
            pilot_stages.main([str(FIXTURES)])
        self.assertEqual(stopped.exception.code, 1)
        self.assertIn("Pilot FID05 refused: sample key is not a pinned pilot key", errors.getvalue())
        self.assertEqual(output.getvalue(), "")
        read_run.assert_not_called()
