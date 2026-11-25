"""Build a local public release repository from a committed Spinivo tree."""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import export_check

NAME = "Erenext LLC"
EMAIL = "spinivo@erenext.com"
MESSAGE = "Initial public release of Spinivo"
GIT_ENV = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
GIT_ENV.update(GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1", GIT_ALLOW_PROTOCOL="",
               GIT_NO_LAZY_FETCH="1", GIT_GRAFT_FILE=os.devnull, GIT_ATTR_NOSYSTEM="1")


def git(root, *args, input=None, stdout=subprocess.PIPE, **env):
    result = subprocess.run(
        ["git", "--no-replace-objects", "-c", "core.hooksPath=/dev/null", "-c", "commit.gpgsign=false",
         "-c", "core.autocrlf=false", *args], cwd=root, env={**GIT_ENV, **env},
        input=input, stdout=stdout, stderr=subprocess.PIPE, check=False)
    if result.returncode:
        raise ValueError(f"git {args[0]} failed: {os.fsdecode(result.stderr).strip()}")
    return result.stdout


def tree_files(root, ref):
    files = {}
    for entry in git(root, "ls-tree", "-r", "-z", "--full-tree", ref).split(b"\0"):
        if not entry:
            continue
        metadata, path = entry.split(b"\t", 1)
        mode, kind, oid = metadata.split()
        name = os.fsdecode(path)
        if mode in (b"120000", b"160000"):
            label = "symlink" if mode == b"120000" else "submodule"
            raise ValueError(f"{label}: {name!r}")
        if kind != b"blob" or mode not in (b"100644", b"100755"):
            raise ValueError(f"unsupported tree entry: {name!r}")
        if name in files or any(part.lower() in ("", ".", "..", ".git") for part in name.split("/")):
            raise ValueError(f"invalid tree path: {name!r}")
        files[name] = mode, oid
    return files


def scan_identities(root, ref, stage, cleared):
    log = git(root, "log", "--no-show-signature", "--no-notes", "--encoding=none", "-z",
              "--format=%an%x00%ae%x00%cn%x00%ce", ref, "--")
    fields = log.lower().removesuffix(b"\0").split(b"\0")
    if not log.endswith(b"\0") or len(fields) % 4:
        raise ValueError("could not read history identities")
    names = set(fields[::2]) - {b"", NAME.lower().encode()}
    emails = set(fields[1::2]) - {b"", EMAIL.lower().encode()}
    identities = names | emails
    if not identities:
        return False
    found = False
    for index, name in enumerate(cleared, 1):
        data = (stage / name).read_bytes().lower()
        for kind, values in (("name", names), ("email", emails)):
            if any(value in data for value in values):
                path = repr(name)
                path_bytes = os.fsencode(name).lower()
                if any(value in path_bytes for value in identities):
                    path = f"cleared file #{index} (path redacted)"
                print(f"history identity {kind}\t{path}", file=sys.stderr)
                found = True
    return found


def check_parent_repo(root, parent_repo, parent_commit, object_format):
    oid_length = {"sha1": 40, "sha256": 64}[object_format]
    if len(parent_commit) != oid_length or not re.fullmatch(r"[0-9a-fA-F]+", parent_commit):
        raise ValueError("--parent-commit must be a full hex commit OID")
    if not parent_repo.is_dir():
        raise ValueError("--parent-repo must be a local repository directory")
    source_dir, parent_dir = (
        Path(os.fsdecode(git(repo, "rev-parse", "--path-format=absolute", "--git-common-dir")
                         .removesuffix(b"\n"))).resolve()
        for repo in (root, parent_repo)
    )
    if source_dir == parent_dir:
        raise ValueError("--parent-repo must be separate from ROOT")
    if git(parent_repo, "rev-parse", "--show-object-format").strip() != object_format.encode():
        raise ValueError("parent and source object formats must match")
    if git(parent_repo, "rev-parse", "--is-shallow-repository").strip() == b"true":
        raise ValueError("--parent-repo needs complete history")


def fetch_parent(output, parent_repo, parent_commit):
    main_oid = git(parent_repo, "show-ref", "--verify", "--hash", "refs/heads/main").strip()
    if main_oid != parent_commit.encode():
        raise ValueError("parent main does not match --parent-commit")
    git(output, "fetch", "--no-tags", "--no-recurse-submodules", "--no-auto-maintenance", "--",
        parent_repo, "refs/heads/main", GIT_ALLOW_PROTOCOL="file")
    actual = git(output, "rev-parse", "--verify", "FETCH_HEAD").strip()
    if actual != parent_commit.encode() or git(output, "cat-file", "-t", actual).strip() != b"commit":
        raise ValueError("parent main does not match --parent-commit")
    return len(verify_history(output, parent_commit))


def commit_release(output, stage, cleared, message, object_format, parent_repo=None, parent_commit=None):
    output.mkdir()
    try:
        git(output, "init", "--template=", "--initial-branch=main", f"--object-format={object_format}")
        parent_count = 0
        if parent_repo is not None:
            parent_count = fetch_parent(output, parent_repo, parent_commit)
        index = []
        for name, (mode, expected) in cleared.items():
            path = output / name
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(stage / name, path)
            path.chmod(0o755 if mode == b"100755" else 0o644)
            oid = git(output, "hash-object", "-w", "--no-filters", "--", path).strip()
            if oid != expected:
                raise ValueError(f"release blob differs from source: {name!r}")
            index.append(mode + b" " + oid + b"\t" + os.fsencode(name) + b"\0")
        git(output, "update-index", "-z", "--index-info", input=b"".join(index))
        date = f"{int(datetime.now(timezone.utc).timestamp())} +0000"
        committer = dict(GIT_COMMITTER_NAME=NAME, GIT_COMMITTER_EMAIL=EMAIL, GIT_COMMITTER_DATE=date)
        if parent_commit is not None:
            if git(output, "write-tree") == git(output, "rev-parse", f"{parent_commit}^{{tree}}"):
                raise ValueError("release tree equals the parent tree")
            git(output, "update-ref", "refs/heads/main", parent_commit, **committer)
        git(output, "commit", "--allow-empty", "--no-gpg-sign", "--cleanup=verbatim", "-m", message,
            GIT_AUTHOR_NAME=NAME, GIT_AUTHOR_EMAIL=EMAIL, GIT_AUTHOR_DATE=date, **committer)
        return verify_release(output, cleared, parent_commit, parent_count)
    except BaseException:
        shutil.rmtree(output)
        raise


def verify_commit(output, commit):
    headers, _, message = git(output, "cat-file", "commit", commit).partition(b"\n\n")
    identity = re.escape(f"{NAME} <{EMAIL}>".encode())
    for role in (b"author", b"committer"):
        fields = re.findall(rb"^" + role + rb" (.*)$", headers, re.MULTILINE)
        if len(fields) != 1 or not re.fullmatch(identity + rb" [0-9]+ \+0000", fields[0]):
            raise ValueError(f"release {role.decode()} identity or UTC offset differs")
    if git(output, "interpret-trailers", "--parse", "--no-divider", input=message).strip():
        raise ValueError("release commit message contains trailers")


def verify_history(output, ref):
    history = [line.split() for line in git(output, "rev-list", "--parents", ref, "--").splitlines()]
    for commit, *parents in history:
        if len(parents) > 1:
            raise ValueError("release history must be linear")
        verify_commit(output, commit)
    return history


def verify_release(output, cleared, parent_commit=None, parent_count=0):
    if int(git(output, "rev-list", "--count", "--all")) != parent_count + 1:
        if parent_commit is None:
            raise ValueError("release must have exactly one commit")
        raise ValueError("release must have exactly one more commit than the parent history")
    history = verify_history(output, "HEAD")
    expected_parents = [] if parent_commit is None else [parent_commit.encode()]
    if history[0][1:] != expected_parents:
        raise ValueError("release HEAD has an unexpected parent")
    if git(output, "symbolic-ref", "HEAD").strip() != b"refs/heads/main":
        raise ValueError("release HEAD must point to main")
    # Same object format as the source, so the IDs verify the exact blob bytes.
    if tree_files(output, "HEAD") != cleared:
        raise ValueError("release tree differs from the cleared paths, modes or blobs")
    commit, tree = git(output, "rev-parse", "HEAD", "HEAD^{tree}").decode().splitlines()
    return commit, tree


def build(root, ref, output, allowlist, message, parent_repo=None, parent_commit=None):
    if (parent_repo is None) != (parent_commit is None):
        raise ValueError("--parent-repo and --parent-commit must be used together")
    if message is None:
        if parent_repo is not None:
            raise ValueError("--message is required for a follow-up release")
        message = MESSAGE
    if os.path.lexists(output):
        raise ValueError("OUTPUT already exists")
    output = output.parent.resolve() / output.name
    if output.is_relative_to(root):
        raise ValueError("OUTPUT must be outside ROOT")
    if Path(os.fsdecode(git(root, "rev-parse", "--show-toplevel").removesuffix(b"\n"))).resolve() != root:
        raise ValueError("ROOT must be the repository root")
    if git(root, "rev-parse", "--is-shallow-repository").strip() == b"true":
        raise ValueError("ROOT needs complete history for the identity scan")
    ref = git(root, "rev-parse", "--verify", "--end-of-options", f"{ref}^{{commit}}").decode().strip()
    object_format = git(root, "rev-parse", "--show-object-format").decode().strip()
    if parent_repo is not None:
        parent_repo = parent_repo.resolve()
        check_parent_repo(root, parent_repo, parent_commit, object_format)
        parent_commit = parent_commit.lower()
    files = tree_files(root, ref)
    with tempfile.TemporaryDirectory() as temp:
        stage = Path(temp)
        for name, (_, oid) in files.items():
            path = stage / name
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("wb") as dest:
                git(root, "cat-file", "blob", oid, stdout=dest)
        policy = stage / "data/redistribution.json"
        if allowlist is not None:
            # Map in-tree paths lexically, so a working-tree symlink cannot swap the committed policy.
            lexical = Path(os.path.abspath(allowlist))
            inside = [path for path in (lexical, lexical.parent.resolve() / lexical.name) if path.is_relative_to(root)]
            policy = stage / inside[0].relative_to(root) if inside else allowlist
        rules, exact_only = export_check.load_allowlist(policy)
        results = export_check.check(stage, files, rules, exact_only)
        failures = [row for row in results if row[0] in export_check.FAILURES]
        cleared = {name: files[name] for status, name, *_ in results if status == "cleared"}
        missing = sorted({"LICENSE", "README.md", "NOTICE"} - cleared.keys())
        for row in failures:
            print("\t".join(row), file=sys.stderr)
        if missing:
            print(f"required release paths not cleared: {', '.join(missing)}", file=sys.stderr)
        if failures or missing:
            return 1
        if scan_identities(root, ref, stage, cleared):
            return 1
        commit, tree = commit_release(output, stage, cleared, message, object_format, parent_repo, parent_commit)
    print(f"cleared {len(cleared)}, excluded {len(files) - len(cleared)}, commit {commit}, tree {tree}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--root", type=Path, default=export_check.PROJECT)
    parser.add_argument("--ref", default="HEAD")
    parser.add_argument("--allowlist", type=Path, help="paths inside ROOT use REF (default: data/redistribution.json)")
    parser.add_argument("--message")
    parser.add_argument("--parent-repo", type=Path, help="local repository with the previous public main")
    parser.add_argument("--parent-commit", help="full OID of the previous public main commit")
    args = parser.parse_args(argv)
    try:
        return build(args.root.resolve(), args.ref, args.output, args.allowlist, args.message,
                     args.parent_repo, args.parent_commit)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
