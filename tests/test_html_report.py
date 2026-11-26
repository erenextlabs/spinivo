import json
import ntpath
import unittest
from html import escape
from html.parser import HTMLParser
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from spinivo import assess_record, load_record
from spinivo.assessment_report import assessment_report, render_text
from spinivo.html_report import render_html
from tests import test_assay_record as records
from tests.test_cli import write_synthetic_bruker


class PageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.stack = []
        self.elements = []
        self.content = []

    def handle_starttag(self, tag, attrs):
        self.elements.append((tag, dict(attrs)))
        if tag != "meta":
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if not self.stack or self.stack[-1] != tag:
            raise ValueError(f"Unbalanced closing tag {tag}: {self.stack}")
        self.stack.pop()

    def handle_data(self, data):
        self.content.append(data)


def report_for(path):
    record = load_record(path)
    return assessment_report(assess_record(record), record)


class HTMLReportTests(unittest.TestCase):
    def assert_page(self, report):
        page = render_html(report)
        parser = PageParser()
        parser.feed(page)
        parser.close()
        self.assertEqual(parser.stack, [])
        self.assertTrue(page.startswith("<!DOCTYPE html>\n"))
        self.assertIn(("meta", {"charset": "utf-8"}), parser.elements)
        self.assertIn('<meta name="viewport"', page)
        self.assertNotIn("<script", page.lower())
        self.assertNotIn("url(", page.lower())
        for _, attrs in parser.elements:
            self.assertTrue({"href", "src", "srcset", "action"}.isdisjoint(attrs))
        content = "".join(parser.content)
        self.assertIn(f'<h1>{escape(report["assay_id"])}</h1>', page)
        self.assertIn(f'Status: {report["status"]}', content)
        for item in (*report["reasons"], *report["not_assessed"]):
            self.assertIn(item, content)
        method = report["method_report"]
        for name in method["inputs"]:
            self.assertIn(f"<td>{escape(name)}</td>", page)
        for row in method["budget"]:
            self.assertIn(f'<td>{escape(row["name"])}</td>', page.split("<h3>Uncertainty budget</h3>")[1])
        self.assertIn(f'{report["schema"]}, spinivo {report["spinivo_version"]}', content)
        self.assertIn("Spinivo reports what the evidence supports. It makes no general accuracy claim.", content)
        checked = next(line.strip() for line in render_text(report).splitlines() if "Checked purity:" in line)
        if report["status"] == "pass":
            self.assertIn(f'<p class="purity">{escape(checked)} (', page)
        else:
            self.assertIn(checked, parser.content)
        return page, parser.content

    def test_fixture_records_withhold_checked_purity(self):
        for name in (records.JEOL, records.EXTERNAL, records.ENTERED):
            with self.subTest(record=name):
                report = report_for(records.RECORDS / f"{name}.json")
                page, content = self.assert_page(report)
                self.assertIn("Checked purity: cannot be assessed", content)
                if name == records.ENTERED:
                    calculation = next(line.strip().split(", u_c = ")[0] for line in render_text(report).splitlines()
                                       if line.strip().startswith("Purity = "))
                    self.assertIn(f"{calculation} (calculation only, not checked)", content)
                    self.assertEqual(page.count(calculation), 1)
                else:
                    self.assertNotIn("Purity = ", page)

    def test_synthetic_pass_shows_checked_purity_and_evidence(self):
        spectral = records.SpectralRecordTests
        spectral.setUpClass()
        data = {**spectral.internal, "assay_id": "synthetic-internal"}
        with TemporaryDirectory() as folder:
            root = Path(folder)
            fid = data["run"]["file"]
            write_synthetic_bruker(root / fid, spectral.raws[fid])
            path = root / "record.json"
            path.write_text(json.dumps(data), encoding="utf-8")
            report = report_for(path)
        self.assertEqual(report["status"], "pass")
        page, content = self.assert_page(report)
        checked = f'Checked purity: {report["purity"]["value"]:.8g} fraction'
        self.assertIn(f"{checked} (99.92 %, U = 0.95 % with k = 2)", content)
        self.assertNotIn("calculation only, not checked", page)
        self.assertIn("<tr><td>relaxation_delay_s</td><td>60</td><td>s</td><td>D[1]</td></tr>", page)
        area = report["runs"]["run"]["regions"][0]["report"]["nominal_area"]["value"]
        self.assertIn(f"<tr><td>area</td><td>pass</td><td>{area:.8g}</td><td>intensity*Hz</td>", page)
        sections = ("Checks", "Baseline shape", "Runs", "Capture pairs", "Method report", "Profiles")
        positions = [page.index(f"<h2>{section}</h2>") for section in sections]
        self.assertEqual(positions, sorted(positions))
        report["method_report"].update(U=None, U_reason="coverage factor is missing", k=None)
        _, content = self.assert_page(report)
        self.assertIn(f"{checked} (99.92 %, U unavailable)", content)

    def test_renders_without_optional_fields(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        del report["runs"]["run"]["run_id"], report["runs"]["run"]["preparation_id"], report["capture_pairs"]
        page, content = self.assert_page(report)
        self.assertIn("Run ID: missing. Preparation ID: missing.", "".join(content))
        self.assertNotIn("<h2>Capture pairs</h2>", page)

    def test_withheld_purity_uses_report_status(self):
        report = report_for(records.RECORDS / f"{records.ENTERED}.json")
        report["status"] = report["method_report"]["status"] = "flag"
        _, content = self.assert_page(report)
        self.assertIn("Checked purity: flag", content)

    def test_escapes_assay_id_and_reason(self):
        report = report_for(records.RECORDS / f"{records.ENTERED}.json")
        payload = '<script>alert("&")</script>'
        report["assay_id"] = payload
        reason = f"Invalid evidence {payload}"
        report["reasons"].append(reason)
        page, content = self.assert_page(report)
        self.assertIn(f"<li>{escape(reason)}</li>", page)
        self.assertIn(reason, content)
        self.assertNotIn(payload, page)

    def test_fixture_run_sources_show_only_last_path_component(self):
        for name in (records.JEOL, records.EXTERNAL):
            with self.subTest(record=name):
                report = report_for(records.RECORDS / f"{name}.json")
                page = render_html(report)
                self.assertNotIn("../", page)
                for run in report["runs"].values():
                    source = Path(run["source"])
                    self.assertIn(f"Source: {escape(source.name)}.", page)
                    self.assertNotIn(escape(str(source.parent)), page)

    def test_private_run_source_is_shortened_in_regions_and_reasons(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        folder = "/fake/private & shared"
        filename = "run<1>.jdf"
        source = f"{folder}/{filename}"
        run = report["runs"]["run"]
        run["source"] = source
        region = run["regions"][0]["report"]
        region["run"] = source
        reason = f"Inspect {source} before sharing."
        region["observations"][0]["reason"] = reason
        report["reasons"].append(reason)
        original = json.dumps(report)

        page = render_html(report)
        self.assertNotIn(escape(folder), page)
        self.assertIn(f"Source: {escape(filename)}.", page)
        self.assertIn(f'<strong>{region["signal"]} {escape(filename)} ', page)
        self.assertIn(f"<li>Inspect {escape(filename)} before sharing.</li>", page)
        self.assertIn(f"<td>Inspect {escape(filename)} before sharing.</td>", page)
        self.assertEqual(json.dumps(report), original)

    def test_windows_run_sources_are_shortened_in_repr_reasons(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        filename = "run.jdf"
        for folder in (r"C:\Users\Private", r"\\server\share"):
            with self.subTest(folder=folder):
                source = f"{folder}\\{filename}"
                report["runs"]["run"]["source"] = source
                report["reasons"] = [f"Inspect {source!r} before sharing."]
                page = render_html(report)
                self.assertIn(f"Source: {filename}.", page)
                reason = f"Inspect {filename!r} before sharing."
                self.assertIn(f"<li>{escape(reason)}</li>", page)
                self.assertNotIn(folder, page)
                self.assertNotIn(repr(folder)[1:-1], page)

    def test_run_source_with_parent_segment_shows_resolved_folder_name(self):
        report = report_for(records.RECORDS / f"{records.EXTERNAL}.json")
        for source in ("exp/10/pdata/..", "10/../10"):
            with self.subTest(source=source):
                report["runs"]["calibrant"]["source"] = source
                report["reasons"] = [f"Inspect {source} before sharing."]
                page = render_html(report)
                self.assertIn("Source: 10.", page)
                self.assertIn("<li>Inspect 10 before sharing.</li>", page)
                self.assertNotIn("exp/10/pdata", page)
                self.assertNotIn(source, page)

    def test_source_without_folder_leaves_other_text_unchanged(self):
        for source in (".", "/", "./run.jdf", "../run.jdf"):
            with self.subTest(source=source):
                report = report_for(records.RECORDS / f"{records.JEOL}.json")
                report["runs"]["run"]["source"] = source
                page = render_html(report)
                self.assertIn(f'spinivo {escape(report["spinivo_version"])}', page)
                self.assertNotIn("1source2", page)

    def test_root_run_source_is_only_shortened_at_token_boundaries(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        report["runs"]["run"]["source"] = "/mL"
        ordinary = r"Units: mg/mL. Other names: /mL.bak, ./mL, //mL, \/mL"
        report["reasons"] = [ordinary, "Inspect /mL."]

        page = render_html(report, hide_folders=())
        self.assertIn("Source: mL. Run ID:", page)
        self.assertIn(f"<li>{ordinary}</li>", page)
        self.assertIn("<li>Inspect mL.</li>", page)

    def test_known_folders_are_hidden_in_run_ids_and_other_text(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        folder = "/home/u/study & records"
        filename = "run<1>.jdf"
        run = report["runs"]["run"]
        run.update(source=f"{folder}/./{filename}", run_id=f"{folder}/{filename}")
        report["assay_id"] = folder
        report["reasons"] = ["/home/user2/x", "/home/u.bak", f"Look in {folder}\\evidence.txt",
                             "Units: mg/mL, version 1.2.0"]
        original = json.dumps(report)

        page = render_html(report, hide_folders=("/", "/home/u", folder))
        self.assertIn(f"Source: {escape(filename)}. Run ID: ./{escape(filename)}.", page)
        self.assertIn("<h1>.</h1>", page)
        self.assertIn("<li>/home/user2/x</li>", page)
        self.assertIn("<li>/home/u.bak</li>", page)
        self.assertIn(r"<li>Look in .\evidence.txt</li>", page)
        self.assertIn("<li>Units: mg/mL, version 1.2.0</li>", page)
        self.assertNotIn(escape(folder), page)
        self.assertEqual(json.dumps(report), original)

    def test_known_folder_is_hidden_before_punctuation(self):
        report = report_for(records.RECORDS / f"{records.ENTERED}.json")
        folder = "/lab/private-experiment"
        report["reasons"] = [f"raw source {folder!r}", f"Inspect {folder}."]

        page = render_html(report, hide_folders=(folder,))
        reason = "raw source '.'"
        self.assertIn(f"<li>{escape(reason)}</li>", page)
        self.assertIn("<li>Inspect ..</li>", page)
        self.assertNotIn(folder, page)

    def test_windows_folder_is_hidden_with_both_slash_spellings(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        folder = r"C:\Users\Ada\study"
        forward = "C:/Users/Ada/study"
        report["runs"]["run"].update(source="run.jdf", run_id=f"{forward}/run.jdf")
        report["reasons"] = [repr(folder), repr(forward)]

        with patch("spinivo.html_report.os", path=ntpath):
            page = render_html(report, hide_folders=(folder,))
        self.assertIn("Source: run.jdf. Run ID: ./run.jdf.", page)
        self.assertEqual(page.count(f"<li>{escape(repr('.'))}</li>"), 2)
        self.assertNotIn(escape(repr(folder)), page)
        self.assertNotIn(forward, page)

    def test_known_folder_used_as_a_run_source_is_hidden(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        folder = "/home/u/private-experiment"
        report["runs"]["run"]["source"] = folder
        report["reasons"] = [f"Inspect {folder}/fid"]
        page = render_html(report, hide_folders=(folder,))
        self.assertIn("Source: .. Run ID:", page)
        self.assertIn("<li>Inspect ./fid</li>", page)

    def test_hides_absolute_resolved_and_repr_folder_spellings(self):
        root = Path(self.enterContext(TemporaryDirectory()))
        target = root / "real folder"
        target.mkdir()
        alias = root / "record's folder"
        alias.symlink_to(target, target_is_directory=True)
        report = report_for(records.RECORDS / f"{records.ENTERED}.json")
        filename = 'pair "evidence".txt'
        report["reasons"] = [str(target / filename), repr(str(alias / filename))]

        page = render_html(report, hide_folders=(alias,))
        self.assertIn(f"<li>{escape(f'./{filename}')}</li>", page)
        self.assertIn(f"<li>{escape(repr(f'./{filename}'))}</li>", page)

    def test_hides_current_and_home_folders_by_default(self):
        report = report_for(records.RECORDS / f"{records.ENTERED}.json")
        report["reasons"] = [f"Work: {Path.cwd()}/record.json", f"Home: {Path.home()}/notes.txt"]
        page = render_html(report)
        marker = "~" if Path.cwd().resolve() == Path.home().resolve() else "."
        self.assertIn(f"<li>Work: {marker}/record.json</li>", page)
        self.assertIn("<li>Home: ~/notes.txt</li>", page)

    def test_short_hidden_folder_keeps_units_and_schema(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        report["reasons"] = ["Molar mass in g/mol", "Inspect /mol/raw.jdf"]
        page = render_html(report, hide_folders=("/mol", "/1"))
        self.assertIn("<li>Molar mass in g/mol</li>", page)
        self.assertIn("<li>Inspect ./raw.jdf</li>", page)
        self.assertIn(escape(report["schema"]), page)

    def test_long_hidden_folder_is_hidden_after_any_prefix(self):
        report = report_for(records.RECORDS / f"{records.JEOL}.json")
        folder = "/home/u/study"
        report["reasons"] = [f"_{folder}/certificate.pdf_", f"-o{folder}/certificate.pdf"]
        page = render_html(report, hide_folders=(folder,))
        self.assertIn("<li>_./certificate.pdf_</li>", page)
        self.assertIn("<li>-o./certificate.pdf</li>", page)
        self.assertNotIn(folder, page)
