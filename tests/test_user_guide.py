import json
import os
import re
import shlex
import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import unquote, urlsplit

from tests import test_assay_record as records
from tests.test_cli import write_synthetic_bruker
from tests.test_readers import PROJECT

GUIDE = PROJECT / "docs" / "user-guide.md"
FENCE = re.compile(r"^```(?P<language>\w+)\n(?P<body>.*?)^```[ \t]*$", re.MULTILINE | re.DOTALL)
EXIT_STATUS = re.compile(r"\(exit status `(\d+)`\)\.\s*$")
EXIT_EXAMPLE = re.compile(r"^\| `(\d+)` \| [^|]+ \| `(spinivo [^`]+)` \|$", re.MULTILINE)
LINK = re.compile(r"\]\(([^)]+)\)")


class UserGuideTests(unittest.TestCase):
    def setUp(self):
        self.guide = GUIDE.read_text(encoding="utf-8")
        self.blocks = list(FENCE.finditer(self.guide))

    def run_cli(self, command):
        name, *args = shlex.split(command)
        self.assertEqual(name, "spinivo")
        return subprocess.run(
            [sys.executable, "-m", "spinivo", *args], cwd=PROJECT, timeout=30,
            env={**os.environ, "PYTHONPATH": "src"}, capture_output=True, text=True,
        )

    def assert_whole_lines(self, excerpt, output):
        self.assertTrue(excerpt.endswith("\n"))
        self.assertIn("\n" + excerpt, "\n" + output)

    def test_cli_examples(self):
        examples = 0
        for index, block in enumerate(self.blocks):
            if block["language"] != "sh":
                continue
            command = block["body"].removesuffix("\n")
            output = self.blocks[index + 1] if index + 1 < len(self.blocks) else None
            has_output = output is not None and output["language"] == "text"
            is_spinivo = command.startswith("spinivo ")
            if not has_output and not is_spinivo:
                continue
            examples += 1
            with self.subTest(command=command):
                self.assertTrue(is_spinivo, "A shell block with text output must use spinivo")
                self.assertTrue(has_output, "A spinivo example must have a text output block")
                self.assertEqual(len(command.splitlines()), 1)
                caption = self.guide[block.end():output.start()]
                status = EXIT_STATUS.search(caption)
                self.assertIsNotNone(status, "The output caption must state the exit status")
                result = self.run_cli(command)
                self.assertEqual(result.returncode, int(status[1]), result.stderr)
                if re.search(r"\bcomplete\b", caption, re.IGNORECASE):
                    self.assertEqual(output["body"], result.stdout)
                else:
                    self.assert_whole_lines(output["body"], result.stdout)
        self.assertGreater(examples, 0)

    def test_exit_status_table(self):
        examples = EXIT_EXAMPLE.findall(self.guide)
        self.assertCountEqual([int(status) for status, _ in examples], [0, 1, 2, 3])
        for status, command in examples:
            with self.subTest(command=command):
                result = self.run_cli(command)
                self.assertEqual(result.returncode, int(status), result.stderr)

    def test_record_excerpts(self):
        record = PROJECT / "tests" / "fixtures" / "records" / "jeol_internal_standard.json"
        excerpts = [block["body"] for block in self.blocks if block["language"] == "json"]
        self.assertGreater(len(excerpts), 0)
        for excerpt in excerpts:
            self.assert_whole_lines(excerpt, record.read_text(encoding="utf-8"))

    def test_synthetic_passing_excerpt(self):
        excerpts = [block["body"] for previous, block in zip(self.blocks, self.blocks[1:])
                    if block["language"] == "text" and previous["language"] != "sh"]
        self.assertEqual(len(excerpts), 1)
        spectral = records.SpectralRecordTests
        spectral.setUpClass()
        data = {**spectral.internal, "assay_id": "synthetic-internal"}
        with TemporaryDirectory() as folder:
            root = Path(folder)
            fid = data["run"]["file"]
            write_synthetic_bruker(root / fid, spectral.raws[fid])
            record = root / "record.json"
            record.write_text(json.dumps(data), encoding="utf-8")
            result = self.run_cli(f"spinivo assess {shlex.quote(str(record))}")
        self.assertEqual((result.returncode, result.stderr), (0, ""))
        self.assert_whole_lines(excerpts[0], result.stdout)

    def test_relative_links(self):
        links = LINK.findall(self.guide)
        self.assertTrue(links)
        for link in links:
            target = urlsplit(link)
            if target.scheme or target.netloc or target.path.startswith("/"):
                continue
            with self.subTest(link=link):
                path = GUIDE.parent / unquote(target.path) if target.path else GUIDE
                self.assertTrue(path.is_file(), path)
        readme = (PROJECT / "README.md").read_text(encoding="utf-8")
        self.assertIn("docs/user-guide.md", LINK.findall(readme))
