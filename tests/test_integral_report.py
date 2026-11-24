import json
import subprocess
import sys
import unittest
from contextlib import redirect_stdout
from dataclasses import asdict, replace
from io import StringIO
from math import inf, nan
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from spinivo.integral_checks import PILOT
from spinivo.integral_components import ComponentEvidence, Line, check_capture_pair, check_region
from spinivo.integral_report import region_report, render_json, render_text
from spinivo.processing import ProcessingSettings, process
from tests.integral_fixtures import (
    FIXTURES, FOREIGN, NOISE_BOUNDS_HZ, SETTINGS, load_cases, make_baseline_raw, make_raw, ppm_bounds, target_lines,
)
from tools.compare_references import PROJECT, main


def fixture_result(name, *, target_source="fixture", thresholds=PILOT, settings=SETTINGS):
    case = next(case for case in load_cases()["cases"] if case["id"] == name)
    bounds, raw = ppm_bounds(case["bounds"]), make_raw(case)
    spectrum = process(raw, settings, regions_ppm=(bounds,))

    def lines(group):
        return tuple(Line(ppm_bounds((hz,))[0], area, width + settings.lb_hz) for hz, area, width in group)

    satellites = case.get("target") == "satellites"
    foreign = FOREIGN.get(case.get("foreign"), ())
    evidence = ComponentEvidence(
        target_lines=lines(target_lines(case)[1:2] if satellites else target_lines(case)), target_source=target_source,
        foreign_lines=lines(foreign), foreign_source="fixture" if foreign else None,
        satellite_convention="include" if satellites else "absent", satellite_source="fixture",
        coupling_hz=130 if satellites else None,
    )
    return check_region(spectrum, bounds, evidence, raw=raw, noise_bounds_ppm=ppm_bounds(NOISE_BOUNDS_HZ),
                        thresholds=thresholds)


class IntegralReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.clean = region_report(fixture_result("clean_singlet"), run="fixture", signal="target", thresholds=PILOT)
        cls.foreign = region_report(fixture_result("foreign_line"), run="fixture", signal="target", thresholds=PILOT)
        cls.missing = region_report(fixture_result("clean_singlet", target_source=None), run="fixture", signal="target",
                                    thresholds=PILOT)

    def test_text_outcomes(self):
        clean, foreign, missing = (render_text(report) for report in (self.clean, self.foreign, self.missing))
        self.assertEqual(clean.splitlines()[0], "target fixture 4.75-5.25 ppm: PASS")
        self.assertIn(f'Area: {self.clean["nominal_area"]["value"]:.8g} intensity*Hz', clean)
        self.assertEqual(foreign.splitlines()[0], "target fixture 4.75-5.25 ppm: FLAG")
        self.assertIn("neighbour: Foreign area fraction exceeds the limit", foreign)
        self.assertEqual(missing.splitlines()[0], "target fixture 4.75-5.25 ppm: CANNOT BE ASSESSED")
        self.assertIn("capture: Missing target source", missing)
        self.assertEqual([report["purity_eligible"] for report in (self.clean, self.foreign, self.missing)],
                         [True, False, False])

    def test_json_observations_and_evidence(self):
        result = fixture_result("satellites_inside")
        report = region_report(result, run="fixture", signal="target", thresholds=PILOT)
        encoded = render_json(report)
        self.assertEqual(encoded, render_json(dict(reversed(list(report.items())))))
        decoded = json.loads(encoded)
        self.assertEqual(decoded, report)
        self.assertEqual(decoded["schema"], "spinivo.integral_region/1")
        self.assertEqual(decoded["thresholds"], asdict(PILOT))
        self.assertEqual(decoded["evidence"], {"target_source": "fixture", "foreign_source": None,
                                               "satellite_source": "fixture", "satellite_convention": "include"})
        self.assertEqual(decoded["nominal_area"], {"value": result.nominal_area, "unit": "intensity*Hz"})
        self.assertEqual(self.foreign["flags"], [
            {"check": "neighbour", "reason": "Foreign area fraction exceeds the limit"},
        ])
        self.assertEqual(self.missing["gaps"], [{"check": "capture", "reason": "Missing target source"}])
        observations = {obs["check"]: obs for obs in decoded["observations"]}
        self.assertEqual(list(observations), [obs.check for obs in result.observations])
        for obs in result.observations:
            with self.subTest(check=obs.check):
                decoded_obs = observations[obs.check]
                self.assertIn("details", decoded_obs)
                for field in ("check", "status", "value", "unit", "limit", "reason", "threshold_version"):
                    self.assertEqual(decoded_obs[field], getattr(obs, field))
        self.assertEqual(observations["bounds"]["details"], json.loads(json.dumps(
            [asdict(trial) for trial in result.bounds.details])))
        self.assertEqual(observations["phase"]["details"], [asdict(trial) for trial in result.phase.details])
        self.assertEqual(observations["satellites"]["details"]["bands"],
                         [asdict(band) for band in result.satellites.bands])
        self.assertEqual(observations["baseline"]["details"], json.loads(json.dumps(asdict(result.baseline.details))))

    def test_spike_anchor_flips_are_reported_without_flagging_the_locked_fit(self):
        bounds = ppm_bounds((-17, 19))
        raw = make_baseline_raw(curvature=0.5)
        settings = ProcessingSettings(zero_fill_factor=1, smooth_baseline=True, baseline_bridge=True)
        unlocked, locked = (process(raw, replace(settings, lock_spike_anchors=lock), regions_ppm=(bounds,))
                            for lock in (False, True))
        self.assertEqual(unlocked.data.tobytes(), locked.data.tobytes())
        self.assertEqual(unlocked.smooth_baseline_curve.tobytes(), locked.smooth_baseline_curve.tobytes())
        before, after = (check_region(spectrum, bounds, ComponentEvidence()) for spectrum in (unlocked, locked))
        self.assertEqual(before.nominal_area.hex(), after.nominal_area.hex())
        self.assertEqual(before.baseline.value, after.baseline.value)
        self.assertEqual((before.baseline.status, after.baseline.status), ("flag", "pass"))
        self.assertGreater(before.baseline.details.recovery.position_sensitivity, 0.017)
        self.assertLess(after.baseline.details.recovery.position_sensitivity, 0.017)
        self.assertNotIn("baseline", {check for check, _ in (*after.flags, *after.gaps)})

        for result, version in ((before, "pilot-2"), (after, "pilot-3")):
            with self.subTest(version=version):
                self.assertEqual({obs.threshold_version for obs in result.observations}, {version})
                report = json.loads(render_json(region_report(result, run="fixture", signal="target", thresholds=PILOT)))
                self.assertEqual(report["thresholds"]["version"], version)
                baseline = next(obs for obs in report["observations"] if obs["check"] == "baseline")
                recovery = baseline["details"]["recovery"]
                if version == "pilot-2":
                    self.assertNotIn("anchor_flips", recovery)
                else:
                    flips = {trial["position_percent"]: trial for trial in recovery["anchor_flips"]}
                    self.assertEqual(set(flips), {25, 50, 75})
                    self.assertEqual(flips[25]["added_ppm"], [5.09125, 4.91125])
                    self.assertEqual(flips[25]["removed_ppm"], [4.93125])
                    self.assertIn("25%: would remove 4.93125 ppm", render_text(report))

    def test_nonfinite_values_keep_reasons(self):
        result = fixture_result("clean_singlet")
        area = replace(result.area, value=nan, status="cannot be assessed", reason="Nominal area is nonfinite")
        trial = replace(result.bounds.details[0], relative_change=inf, reason="Relative area change is nonfinite")
        result = replace(result, nominal_area=nan, area=area, status=area.status, gaps=(("area", area.reason),),
                         bounds=replace(result.bounds, details=(trial, *result.bounds.details[1:])))
        decoded = json.loads(render_json(region_report(result, run="fixture", signal="target", thresholds=PILOT)))
        observations = {obs["check"]: obs for obs in decoded["observations"]}
        self.assertIsNone(decoded["nominal_area"]["value"])
        self.assertIsNone(observations["area"]["value"])
        self.assertEqual(observations["area"]["reason"], area.reason)
        self.assertIsNone(observations["bounds"]["details"][0]["relative_change"])
        self.assertEqual(observations["bounds"]["details"][0]["reason"], trial.reason)

    def test_custom_threshold_profile(self):
        thresholds = replace(PILOT, min_noise_points=64)
        result = fixture_result("clean_singlet", thresholds=thresholds)
        report = region_report(result, run="fixture", signal="target", thresholds=thresholds)
        self.assertEqual(report["thresholds"], asdict(thresholds))
        for change in ({"version": "test"}, {"max_relative_change": 0.02}, {"min_snr": 200}):
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, "does not match"):
                region_report(result, run="fixture", signal="target", thresholds=replace(thresholds, **change))
        for check, change in (("capture", {"limit": 0.02}), ("neighbour", {"threshold_version": "other"}),
                              ("phase", {"limit": 0.02}), ("snr", {"threshold_version": "other"})):
            recorded = replace(result, **{check: replace(getattr(result, check), **change)})
            with self.subTest(check=check), self.assertRaisesRegex(ValueError, "does not match"):
                region_report(recorded, run="fixture", signal="target", thresholds=thresholds)

    def test_profile_compares_limits_without_observation_fields(self):
        result = fixture_result("clean_singlet")
        for change in ({"unexpected_peak_sd": 6}, {"min_noise_points": 33}, {"max_pair_capture_mismatch": 0.006}):
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, "does not match"):
                region_report(result, run="fixture", signal="target", thresholds=replace(PILOT, **change))

    def test_capture_pair_requires_matching_applied_anchor_locks(self):
        settings = replace(SETTINGS, lock_spike_anchors=True)
        locked = fixture_result("clean_singlet", settings=settings)
        pair = check_capture_pair(locked, locked)
        self.assertEqual((pair.status, pair.threshold_version), ("pass", "pilot-3"))

        windows = (ppm_bounds((-300, -200)), ppm_bounds((200, 300)))
        linear = replace(settings, smooth_baseline=False, baseline_windows_ppm=windows)
        unlocked = fixture_result("clean_singlet", settings=linear)
        self.assertEqual({obs.threshold_version for obs in unlocked.observations}, {"pilot-2"})
        pair = check_capture_pair(unlocked, unlocked)
        self.assertEqual((pair.status, pair.threshold_version), ("pass", "pilot-2"))
        for analyte, standard in ((locked, unlocked), (unlocked, locked)):
            pair = check_capture_pair(analyte, standard)
            self.assertEqual((pair.status, pair.value), ("cannot be assessed", None))
            self.assertIn("different spike anchor lock modes", pair.reason)

    def test_locked_region_preserves_custom_threshold_profile(self):
        thresholds = replace(PILOT, version="fixture", min_noise_points=64)
        settings = replace(SETTINGS, lock_spike_anchors=True)
        result = fixture_result("clean_singlet", thresholds=thresholds, settings=settings)
        self.assertEqual({obs.threshold_version for obs in result.observations}, {thresholds.version})
        report = region_report(result, run="fixture", signal="target", thresholds=thresholds)
        self.assertEqual(report["thresholds"], asdict(thresholds))
        self.assertEqual(check_capture_pair(result, result, thresholds).threshold_version, thresholds.version)


class CompareChecksTests(unittest.TestCase):
    def test_checks_cli(self):
        completed = subprocess.run(
            [sys.executable, str(PROJECT / "tools/compare_references.py"), "--checks", str(FIXTURES)],
            cwd=FIXTURES, check=True, capture_output=True, text=True,
        )
        lines = completed.stdout.splitlines()
        for name, bounds, label in (("bruker_1_15", "1.15-1.22", "unreliable"),
                                    ("bruker_7_64", "7.64-7.69", "unreliable"),
                                    ("fid02_dmso2", "2.9-3.01", "indeterminate"),
                                    ("fid02_ma", "6.17-6.28", "indeterminate")):
            with self.subTest(case=name):
                start = next(i for i, line in enumerate(lines) if line.startswith(f"{name} "))
                end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("label="))
                self.assertEqual(lines[end], f"label={label}")
                self.assertIn(f" {bounds} ppm: ", lines[start])
                status = lines[start].rsplit(": ", 1)[1]
                body = "\n".join(lines[start + 1:end])
                self.assertIn("capture: Missing target source", body)
                if name.startswith("bruker"):
                    self.assertEqual(status, "FLAG")
                    self.assertRegex(body, r"bounds: Relative area change exceeds the limit\n"
                                           r"(?:    (?:low|high)_(?:in|out) [^\n]+ ppm: [+-][\d.]+%\n){4}")
                else:
                    self.assertIn(status, ("FLAG", "CANNOT BE ASSESSED"))
        self.assertEqual(sum(line.startswith("label=") for line in lines), 4)

    def test_checks_cli_deposit_path_and_shared_run(self):
        cases = {case["id"]: case for case in load_cases()["cases"] if case["kind"] == "real"}
        widened = cases["fid05_widened_dmso2"]
        output = StringIO()
        with TemporaryDirectory() as root:
            path = Path(root) / widened["source"]
            path.parent.mkdir(parents=True)
            path.symlink_to(PROJECT / cases["fid02_dmso2"]["source"])
            with (patch.object(sys, "argv", ["compare_references.py", "--checks", root]),
                  patch("spinivo.processing.process", wraps=process) as process_run,
                  redirect_stdout(output)):
                main()

        lines = output.getvalue().splitlines()
        headings = [line for line in lines if " ppm: " in line and not line.startswith(" ")]
        self.assertEqual([line.split()[0] for line in headings], list(cases))
        self.assertTrue(headings[-1].startswith(f'fid05_widened_dmso2 {widened["source"]} 2.93-3.48 ppm: '))
        self.assertEqual(lines[-1], "label=unreliable")
        self.assertEqual([call.kwargs["regions_ppm"] for call in process_run.call_args_list], [
            ((2.90, 3.01), (6.17, 6.28), (2.93, 3.48)), ((1.15, 1.22), (7.64, 7.69)),
        ])
