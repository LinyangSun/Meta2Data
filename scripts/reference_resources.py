#!/usr/bin/env python3
"""Prepare shared QIIME references in <working directory>/results/db.

This internal helper is shared by AmpliconPIP and AmpliconTAXA. Its CLI prints
only the prepared artifact path to stdout; progress and errors go to stderr.
"""
import argparse
from dataclasses import dataclass
import fcntl
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import uuid
import zipfile
import zlib


@dataclass(frozen=True)
class Resource:
    filename: str
    semantic_type: str
    url: str


GG2_BASE = "https://ftp.microbio.me/greengenes_release/2024.09/"
RESOURCES = {
    "gg2-sequences": Resource(
        "2024.09.backbone.full-length.fna.qza", "FeatureData[Sequence]",
        GG2_BASE + "2024.09.backbone.full-length.fna.qza"),
    "greengenes-classifier": Resource(
        "2024.09.backbone.full-length.nb.qza", "TaxonomicClassifier",
        GG2_BASE + "2024.09.backbone.full-length.nb.qza"),
    "silva-classifier": Resource(
        "silva-138-99-nb-classifier.qza", "TaxonomicClassifier",
        "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"),
    "sepp": Resource(
        "sepp-refs-gg-13-8.qza", "SeppReferenceDatabase",
        "https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-gg-13-8.qza"),
}


def database_directory(work_dir=None):
    """Resolve the shared cache independently of each module's output directory."""
    return Path(work_dir or Path.cwd()).resolve() / "results" / "db"


def validate_artifact(path, semantic_type):
    """Check archive integrity, UUID, semantic type, and a nonempty data payload."""
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        metadata_paths = [name for name in names
                          if name.count("/") == 1 and name.endswith("/metadata.yaml")]
        if len(metadata_paths) != 1:
            raise ValueError("expected one root metadata.yaml")
        metadata = {}
        for line in archive.read(metadata_paths[0]).decode("utf-8").splitlines():
            if ":" in line and not line[:1].isspace():
                key, value = line.split(":", 1)
                metadata[key] = value.strip().strip("'\"")
        artifact_uuid = str(uuid.UUID(metadata.get("uuid", "")))
        if artifact_uuid != metadata_paths[0].split("/")[0]:
            raise ValueError("artifact UUID does not match archive root")
        if metadata.get("type") != semantic_type:
            raise ValueError(f"expected {semantic_type}, found {metadata.get('type')}")
        if not any(info.filename.startswith(f"{artifact_uuid}/data/")
                   and not info.is_dir() and info.file_size > 0 for info in archive.infolist()):
            raise ValueError("artifact contains no nonempty data files")
        bad_member = archive.testzip()
        if bad_member:
            raise ValueError(f"archive integrity check failed: {bad_member}")


def valid_artifact(path, semantic_type):
    try:
        validate_artifact(path, semantic_type)
        return True
    except (OSError, ValueError, UnicodeError, zipfile.BadZipFile, RuntimeError, EOFError, zlib.error):
        return False


def bundled_candidates(kind, package_root=None, environ=None):
    """Find only explicit deployment resources and conventional package locations."""
    environ = os.environ if environ is None else environ
    package_root = Path(package_root or Path(__file__).resolve().parent.parent)
    if kind == "sepp" and environ.get("META2DATA_SEPP_REF"):
        yield Path(environ["META2DATA_SEPP_REF"])
    roots = [package_root,
             Path("/usr/local/share/Meta2Data"),
             Path("/opt/Meta2Data"), Path("/opt/meta2data")]
    for variable in ("CONDA_PREFIX", "PREFIX"):
        if environ.get(variable):
            roots.append(Path(environ[variable]) / "share" / "Meta2Data")
    seen = set()
    for root in roots:
        for subdirectory in ("references", "reference", "refs", "db"):
            candidate = root / subdirectory / RESOURCES[kind].filename
            if str(candidate) not in seen:
                seen.add(str(candidate))
                yield candidate


def download(url, destination):
    """Download with certificate verification; never publish a partial artifact."""
    subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error",
                    "--retry", "3", "--connect-timeout", "30",
                    "--output", str(destination), url], check=True)


def prepare_resource(kind, work_dir=None, *, builtin_paths=None, downloader=None):
    """Return a valid cached artifact, copying a bundle or downloading if needed.

    builtin_paths and downloader are internal integration/test hooks.
    All callers cache resources beneath the startup working directory.
    """
    if kind not in RESOURCES:
        raise ValueError(f"Unknown reference resource: {kind}")
    resource = RESOURCES[kind]
    directory = database_directory(work_dir)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / resource.filename
    with (directory / f".{resource.filename}.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if valid_artifact(target, resource.semantic_type):
            return target
        if target.exists():
            print(f"Replacing invalid reference: {target}", file=sys.stderr)
        source = None
        candidates = bundled_candidates(kind) if builtin_paths is None else builtin_paths
        for candidate in candidates:
            candidate = Path(candidate)
            if candidate.resolve() == target.resolve():
                continue
            if valid_artifact(candidate, resource.semantic_type):
                source = candidate
                break
            if candidate.is_file():
                print(f"Ignoring invalid bundled reference: {candidate}", file=sys.stderr)
        with tempfile.NamedTemporaryFile(prefix=f".{resource.filename}.", suffix=".partial",
                                         dir=directory, delete=False) as temporary:
            staging = Path(temporary.name)
        try:
            if source is not None:
                print(f"Preparing bundled reference: {source} -> {target}", file=sys.stderr)
                shutil.copyfile(source, staging)
            else:
                print(f"Downloading {kind} to {target}", file=sys.stderr)
                (downloader or download)(resource.url, staging)
            validate_artifact(staging, resource.semantic_type)
            os.replace(staging, target)
        except (OSError, ValueError, UnicodeError, zipfile.BadZipFile, RuntimeError,
                EOFError, zlib.error, subprocess.CalledProcessError) as error:
            raise RuntimeError(
                f"Could not prepare {kind} in {directory}: {error}. "
                f"Check network access, free disk space, and write permission; rerun to retry. "
                f"Reference URL: {resource.url}") from error
        finally:
            staging.unlink(missing_ok=True)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare", help="Reuse, copy, or download a valid reference")
    prepare.add_argument("--kind", choices=RESOURCES, required=True)
    prepare.add_argument("--work-dir", default=None, help="Module startup directory (default: cwd)")
    args = parser.parse_args()
    try:
        print(prepare_resource(args.kind, work_dir=args.work_dir))
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
