import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

TOOL = Path(__file__).resolve().parents[1] / "tools" / "build_release_tree.py"
IDENTITY = b"Erenext LLC <spinivo@erenext.com>"
FILE = "src/run\tfixture.bin"
CONTENTS = b"#!/bin/sh\r\n\xff\x00\r\n"
AUTHOR = "Synthetic Author"
COMMITTER_EMAIL = "fixture-committer@example.invalid"
REQUIRED = ("LICENSE", "NOTICE", "README.md")


class BuildReleaseTreeTests(unittest.TestCase):
    def setUp(self):
        self.temp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.root = self.temp / "source"
        self.root.mkdir()
        self.output = self.temp / "release"
        self.allowlist = self.temp / "redistribution.json"
        home = self.temp / "home"
        home.mkdir()
        self.env = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
        self.env.update(HOME=str(home), XDG_CONFIG_HOME=str(home / ".config"), TZ="Etc/GMT-9")
        self.git_env = dict(self.env, GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1",
                            GIT_AUTHOR_NAME=AUTHOR, GIT_AUTHOR_EMAIL="fixture-author@example.invalid",
                            GIT_COMMITTER_NAME="Synthetic Committer", GIT_COMMITTER_EMAIL=COMMITTER_EMAIL)
        self.git("init", "--template=", "--initial-branch=main")
        self.write(FILE, CONTENTS).chmod(0o755)
        self.write("private.txt", AUTHOR.encode())
        for name in REQUIRED:
            self.write(name, b"Synthetic release documentation\n")
        self.git("add", "--all")
        self.git("commit", "-m", "Synthetic source")
        self.rules = {
            **{name: self.rule(name) for name in REQUIRED},
            FILE: self.rule(FILE, sha256=hashlib.sha256(CONTENTS).hexdigest()),
            "private.txt": self.rule("private.txt", distributable=False),
        }

    def git(self, *args, root=None, input=None, env=None):
        result = subprocess.run(["git", "-c", "core.filemode=true", "-c", "core.autocrlf=false", *args],
                                cwd=root or self.root, input=input, env=env or self.git_env,
                                capture_output=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        return result.stdout

    def write(self, name, data):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def rule(self, path, **changes):
        return {"path": path, "licence": "CC0-1.0", "source": "synthetic", "modifications": "none",
                "reason": "test fixture", "distributable": True, **changes}

    def run_builder(self, *args, output=None, use_external_allowlist=True):
        self.allowlist.write_text(json.dumps({"rules": list(self.rules.values())}), encoding="utf-8")
        policy = ["--allowlist", self.allowlist] if use_external_allowlist else []
        return subprocess.run([sys.executable, TOOL, output or self.output, "--root", self.root, *policy, *args],
                              cwd=self.temp, env=self.env, capture_output=True, check=False)

    def commit_policy(self):
        name = "data/redistribution.json"
        self.rules[name] = self.rule(name)
        path = self.write(name, json.dumps({"rules": list(self.rules.values())}).encode())
        self.git("add", name)
        self.git("commit", "-m", "Record policy")
        return path

    def assert_failed(self, result, details, output=None):
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(details, result.stderr)
        self.assertFalse(os.path.lexists(output or self.output))

    def test_release_ignores_hostile_configuration_and_attributes(self):
        hook = self.temp / "hooks" / "commit-msg"
        hook.parent.mkdir()
        hook.write_text("#!/bin/sh\ntouch ../hook-ran\nprintf '\\nCo-authored-by: Hook User\\n' >> \"$1\"\n",
                        encoding="utf-8")
        hook.chmod(0o755)
        (Path(self.env["HOME"]) / ".gitconfig").write_text(
            f"[user]\n\tname = Wrong User\n\temail = wrong@example.invalid\n"
            f"[commit]\n\tgpgsign = true\n[gpg]\n\tprogram = /no-such-gpg\n[core]\n\thooksPath = {hook.parent}\n",
            encoding="utf-8")
        self.env.update(GIT_CONFIG_COUNT="1", GIT_CONFIG_KEY_0="commit.gpgsign", GIT_CONFIG_VALUE_0="true",
                        GIT_DIR=str(self.root / ".git"), GIT_AUTHOR_DATE="1000000000 +0900")
        self.write(".gitattributes", b"# " + IDENTITY + b"\n*.bin text\n")
        self.git("add", ".gitattributes")
        release_env = dict(self.git_env, GIT_AUTHOR_NAME="Erenext LLC", GIT_AUTHOR_EMAIL="spinivo@erenext.com",
                           GIT_COMMITTER_NAME="Erenext LLC", GIT_COMMITTER_EMAIL="spinivo@erenext.com")
        self.git("commit", "-m", "Add attributes after recording raw bytes", env=release_env)
        self.rules[".gitattributes"] = self.rule(".gitattributes")

        result = self.run_builder()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        commit = self.git("rev-parse", "HEAD", root=self.output).strip()
        tree = self.git("rev-parse", "HEAD^{tree}", root=self.output).strip()
        self.assertEqual(result.stdout, b"cleared 5, excluded 1, commit %s, tree %s\n" % (commit, tree))
        self.assertEqual(self.git("rev-list", "--count", "--all", root=self.output), b"1\n")
        self.assertEqual(self.git("symbolic-ref", "HEAD", root=self.output), b"refs/heads/main\n")
        self.assertEqual(self.git("remote", root=self.output), b"")
        raw = self.git("cat-file", "commit", "HEAD", root=self.output)
        headers, _, message = raw.partition(b"\n\n")
        for role in (b"author", b"committer"):
            self.assertRegex(headers, rb"(?m)^" + role + rb" " + IDENTITY + rb" [0-9]+ \+0000$")
        self.assertNotIn(b"gpgsig", headers)
        self.assertEqual(message, b"Initial public release of Spinivo\n")
        self.assertEqual(self.git("interpret-trailers", "--parse", root=self.output, input=message), b"")
        self.assertFalse((self.temp / "hook-ran").exists())
        self.assertFalse((self.output / ".git" / "hooks").exists())
        self.assertFalse((self.output / "private.txt").exists())
        self.assertEqual(self.git("ls-tree", "-r", "-z", "--name-only", "HEAD", root=self.output),
                         b".gitattributes\0LICENSE\0NOTICE\0README.md\0" + os.fsencode(FILE) + b"\0")
        self.assertEqual(self.git("ls-tree", "-r", "-z", "HEAD", "--", FILE, root=self.output),
                         self.git("ls-tree", "-r", "-z", "HEAD", "--", FILE))
        self.assertEqual(self.git("cat-file", "blob", f"HEAD:{FILE}", root=self.output), CONTENTS)
        self.assertEqual((self.output / FILE).read_bytes(), CONTENTS)
        self.assertEqual((self.output / FILE).stat().st_mode & 0o111, 0o111)

    def test_unknown_mismatch_and_missing_credit(self):
        original = self.rules[FILE]
        for status, rule in ((b"unknown", None), (b"mismatch", dict(original, sha256="0" * 64)),
                             (b"attribution", dict(original, licence="CC-BY-4.0", credit_files=["CREDIT"]))):
            with self.subTest(status=status):
                self.rules.pop(FILE, None)
                if rule is not None:
                    self.rules[FILE] = rule
                self.assert_failed(self.run_builder(), status + b"\t")

    def test_ref_without_license(self):
        self.git("rm", "--cached", "LICENSE")
        self.git("commit", "-m", "Remove licence from REF")
        self.assert_failed(self.run_builder(), b"required release paths not cleared: LICENSE")

    def test_symlink_and_submodule(self):
        (self.root / "link").symlink_to("private.txt")
        self.git("add", "link")
        self.git("commit", "-m", "Add symlink")
        self.assert_failed(self.run_builder(), b"symlink")
        self.git("rm", "link")
        oid = self.git("rev-parse", "HEAD").strip().decode()
        self.git("update-index", "--add", "--cacheinfo", f"160000,{oid},submodule")
        self.git("commit", "-m", "Add submodule entry")
        self.assert_failed(self.run_builder(), b"submodule")

    def test_output_location(self):
        self.output.mkdir()
        keep = self.output / "keep"
        keep.write_bytes(b"existing output")
        result = self.run_builder()
        self.assertEqual(result.returncode, 1)
        self.assertIn(b"OUTPUT already exists", result.stderr)
        self.assertEqual(keep.read_bytes(), b"existing output")
        dangling = self.temp / "dangling"
        target = self.temp / "absent"
        dangling.symlink_to(target)
        result = self.run_builder(output=dangling)
        self.assertEqual(result.returncode, 1)
        self.assertIn(b"OUTPUT already exists", result.stderr)
        self.assertFalse(target.exists())
        alias = self.temp / "alias"
        alias.symlink_to(self.root, target_is_directory=True)
        for inside in (self.root / "release", alias / "release"):
            self.assert_failed(self.run_builder(output=inside), b"OUTPUT must be outside ROOT", output=inside)

    def test_ref_ignores_uncommitted_files_and_policy(self):
        policy = self.commit_policy()
        self.git("tag", "candidate")
        self.write(FILE, b"new committed bytes")
        self.git("commit", "-am", "Advance HEAD")
        self.write(FILE, b"working file")
        policy.write_bytes(b"working policy")
        self.write("untracked", b"working file")
        result = self.run_builder("--ref", "candidate", "--message", "Candidate", use_external_allowlist=False)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual((self.output / FILE).read_bytes(), CONTENTS)
        self.assertFalse((self.output / "untracked").exists())
        self.assertEqual(self.git("log", "-1", "--format=%s", root=self.output), b"Candidate\n")

    def test_explicit_in_tree_allowlist_uses_ref(self):
        policy = self.commit_policy()
        self.rules[FILE]["distributable"] = False
        policy.write_text(json.dumps({"rules": list(self.rules.values())}), encoding="utf-8")
        explicit = self.temp / "explicit"
        for root, args in ((self.output, ()), (explicit, ("--allowlist", policy))):
            result = self.run_builder(*args, output=root, use_external_allowlist=False)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.git("rev-parse", "HEAD^{tree}", root=explicit),
                         self.git("rev-parse", "HEAD^{tree}", root=self.output))

    def test_explicit_in_tree_allowlist_symlink_uses_ref(self):
        policy = self.commit_policy()
        self.rules["private.txt"]["distributable"] = True
        swapped = self.temp / "swapped.json"
        swapped.write_text(json.dumps({"rules": list(self.rules.values())}), encoding="utf-8")
        policy.unlink()
        policy.symlink_to(swapped)
        result = self.run_builder("--allowlist", policy, use_external_allowlist=False)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse((self.output / "private.txt").exists())

    def test_history_identity_is_redacted(self):
        self.rules[FILE] = self.rule(FILE)
        self.write(FILE, b"\xff" + AUTHOR.swapcase().encode() + b"\x00" + COMMITTER_EMAIL.upper().encode())
        env = dict(self.git_env, GIT_AUTHOR_NAME="New Author", GIT_AUTHOR_EMAIL="new-author@example.invalid",
                   GIT_COMMITTER_NAME="New Committer", GIT_COMMITTER_EMAIL="new-committer@example.invalid")
        self.git("commit", "-am", "Add bytes containing older identities", env=env)
        result = self.run_builder()
        self.assert_failed(result, b"history identity name\t" + repr(FILE).encode())
        self.assertIn(b"history identity email\t" + repr(FILE).encode(), result.stderr)
        for identity in (AUTHOR, COMMITTER_EMAIL):
            self.assertNotIn(identity.lower().encode(), (result.stdout + result.stderr).lower())

    def test_message_trailers_remove_output(self):
        result = self.run_builder("--message", "Release\n\nSigned-off-by: Synthetic <signer@example.invalid>")
        self.assert_failed(result, b"release commit message contains trailers")
