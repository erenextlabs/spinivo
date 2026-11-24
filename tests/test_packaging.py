import ast
import os
import subprocess
import sys
import tomllib
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO

from spinivo import __version__
from spinivo.cli import main
from tests.test_readers import JEOL, PROJECT


class PackagingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (PROJECT / "pyproject.toml").open("rb") as source:
            cls.config = tomllib.load(source)
        cls.project = cls.config["project"]

    def test_version_comes_from_one_string_literal(self):
        self.assertNotIn("version", self.project)
        self.assertEqual(self.project["dynamic"], ["version"])
        self.assertEqual(self.config["tool"]["setuptools"]["dynamic"]["version"],
                         {"attr": "spinivo.__version__"})

        path = PROJECT / "src" / "spinivo" / "__init__.py"
        module = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        targets = [node for node in ast.walk(module)
                   if isinstance(node, ast.Name) and node.id == "__version__" and isinstance(node.ctx, ast.Store)]
        self.assertEqual(len(targets), 1)
        assignments = [node for node in module.body
                       if isinstance(node, ast.Assign) and node.targets == targets]
        self.assertEqual(len(assignments), 1)
        value = assignments[0].value
        self.assertIsInstance(value, ast.Constant)
        self.assertIsInstance(value.value, str)
        self.assertRegex(value.value, r"\A\d+\.\d+\.\d+((a|b|rc)\d+)?(\.dev\d+)?\Z")
        self.assertEqual(value.value, __version__)

    def test_project_metadata(self):
        self.assertEqual(self.project["name"], "spinivo")
        self.assertEqual(self.project["authors"], [{"name": "Erenext LLC", "email": "spinivo@erenext.com"}])
        repo = "https://github.com/erenextlabs/spinivo"
        self.assertEqual(self.project["urls"], {"Homepage": repo, "Source": repo, "Issues": f"{repo}/issues"})
        self.assertEqual(self.project["requires-python"], ">=3.11")
        self.assertIn("nmrglue>=0.12", self.project["dependencies"])
        self.assertIn("numpy>=1.23.5", self.project["dependencies"])
        classifiers = self.project["classifiers"]
        self.assertIn("Programming Language :: Python :: 3 :: Only", classifiers)
        self.assertIn("Programming Language :: Python :: 3.11", classifiers)
        self.assertFalse(any(item.startswith("License ::") for item in classifiers))

    def test_distribution_documents(self):
        self.assertEqual(self.project["readme"], "README.md")
        self.assertEqual(self.project["license"], "GPL-3.0-or-later")
        self.assertEqual(self.project["license-files"], ["LICENSE", "NOTICE"])
        for name in ("README.md", "LICENSE", "NOTICE"):
            with self.subTest(name=name):
                self.assertTrue((PROJECT / name).read_text(encoding="utf-8").strip())
        notice = (PROJECT / "NOTICE").read_text(encoding="utf-8")
        for detail in ("Spinivo", "Copyright 2026 Erenext LLC", "GPL-3.0-or-later",
                       "Jingyi Hu", "qNMRrawdata_EnhancingPolyketideProduction",
                       "10.5281/zenodo.17083860", "CC BY 4.0",
                       "Nishizaki et al.", "10.7910/DVN/QFY6L9", "CC0 1.0"):
            with self.subTest(detail=detail):
                self.assertIn(detail, notice)

    def test_module_matches_console_entry_point(self):
        self.assertEqual(self.project["scripts"]["spinivo"], "spinivo.cli:main")
        args = ["inspect", str(JEOL), "--format", "json"]
        stdout, stderr = StringIO(), StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            status = main(args)
        self.assertEqual((status, stderr.getvalue()), (0, ""))
        self.assertTrue(stdout.getvalue())

        result = subprocess.run(
            [sys.executable, "-m", "spinivo", *args], cwd=PROJECT, timeout=30,
            env={**os.environ, "PYTHONPATH": "src"}, capture_output=True, text=True,
        )
        self.assertEqual((result.returncode, result.stdout, result.stderr), (status, stdout.getvalue(), ""))
