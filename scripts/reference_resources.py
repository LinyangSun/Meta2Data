#!/usr/bin/env python3
"""Prepare shared QIIME references in <working directory>/results/db.

This internal helper is shared by AmpliconPIP and AmpliconTAXA. Its CLI prints
only the prepared artifact path to stdout; progress and errors go to stderr.
"""
import argparse
from dataclasses import dataclass
import fcntl
import hashlib
import json
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
    url: str | None
    classifier_id: str | None = None
    sklearn_version: str | None = None
    expected_uuid: str | None = None
    md5: str | None = None
    sha256: str | None = None


GG2_BASE = "https://ftp.microbio.me/greengenes_release/2024.09/"
CLASSIFIER_SKLEARN_VERSION = "1.7.1"
SILVA144_BASE = (
    "https://www.arb-silva.de/fileadmin/silva_databases/release_144/"
    "QIIME2/2026.7/SSU/full-length/uniform/")
RESOURCES = {
    "gg2-sequences": Resource(
        "2024.09.backbone.full-length.fna.qza", "FeatureData[Sequence]",
        GG2_BASE + "2024.09.backbone.full-length.fna.qza"),
    "greengenes-classifier": Resource(
        "gg2-2024.09-full-length-sklearn-1.7.1.qza", "TaxonomicClassifier",
        None, classifier_id="gg2-2024.09", sklearn_version=CLASSIFIER_SKLEARN_VERSION),
    "silva-classifier": Resource(
        "silva-138.1-full-length-sklearn-1.7.1.qza", "TaxonomicClassifier",
        None, classifier_id="silva-138.1", sklearn_version=CLASSIFIER_SKLEARN_VERSION),
    "silva144-classifier": Resource(
        "silva-144-full-length-sklearn-1.7.1.qza", "TaxonomicClassifier",
        SILVA144_BASE +
        "SILVA_144_SSURef_NR99_full-length_prokaryotes_rank-propagation_uniq_genus_uniform_classifier.qza",
        classifier_id="silva-144", sklearn_version=CLASSIFIER_SKLEARN_VERSION,
        expected_uuid="b7986601-1851-4fa2-ace5-98bacf660fa3",
        md5="e16b9bc78dafa2f344563e570bbf25ba"),
    "sepp": Resource(
        "sepp-refs-gg-13-8.qza", "SeppReferenceDatabase",
        "https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-gg-13-8.qza"),
}

# The two legacy classifiers must be rebuilt with 1.7.1 before this development
# branch is released. A missing URL is intentional: never substitute the old
# 1.4.2 model or a different SILVA release. Maintainer-built bundles can be used
# while their final UUIDs, checksums, and published URLs are being established.


def database_directory(work_dir=None):
    """Resolve the shared cache independently of each module's output directory."""
    return Path(work_dir or Path.cwd()).resolve() / "results" / "db"


def validate_artifact(path, semantic_type, *, sklearn_version=None, expected_uuid=None,
                      require_nonempty_data=True, check_crc=True):
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
        if expected_uuid is not None and artifact_uuid != expected_uuid:
            raise ValueError(f"expected artifact UUID {expected_uuid}, found {artifact_uuid}")
        if metadata.get("type") != semantic_type:
            raise ValueError(f"expected {semantic_type}, found {metadata.get('type')}")
        data_files = [info for info in archive.infolist()
                      if info.filename.startswith(f"{artifact_uuid}/data/") and not info.is_dir()]
        if not data_files or (require_nonempty_data and not any(info.file_size for info in data_files)):
            raise ValueError("artifact contains no nonempty data files")
        if semantic_type == "TaxonomicClassifier":
            version_path = f"{artifact_uuid}/data/sklearn_version.json"
            try:
                stored_version = json.loads(archive.read(version_path))["sklearn-version"]
            except (KeyError, TypeError, json.JSONDecodeError) as error:
                raise ValueError("classifier is missing valid sklearn_version.json metadata") from error
            if not isinstance(stored_version, str) or not stored_version:
                raise ValueError("classifier has no recorded scikit-learn version")
            if sklearn_version is not None and stored_version != sklearn_version:
                raise ValueError(
                    f"classifier requires scikit-learn {stored_version}; "
                    f"this Meta2Data build requires {sklearn_version}")
            pipeline = f"{artifact_uuid}/data/sklearn_pipeline.tar"
            if pipeline not in names or archive.getinfo(pipeline).file_size == 0:
                raise ValueError("classifier contains no trained model payload")
        if check_crc:
            bad_member = archive.testzip()
            if bad_member:
                raise ValueError(f"archive integrity check failed: {bad_member}")


def valid_artifact(path, semantic_type, **kwargs):
    try:
        validate_artifact(path, semantic_type, **kwargs)
        return True
    except (OSError, ValueError, UnicodeError, zipfile.BadZipFile, RuntimeError, EOFError, zlib.error):
        return False


def file_digest(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_resource(path, resource):
    """Validate a fixed built-in model without loading or unpickling it."""
    # Pinned digests below replace the redundant full-payload CRC pass.
    validate_artifact(path, resource.semantic_type,
                      sklearn_version=resource.sklearn_version,
                      expected_uuid=resource.expected_uuid,
                      check_crc=resource.md5 is None and resource.sha256 is None)
    for algorithm in ("md5", "sha256"):
        expected = getattr(resource, algorithm)
        if expected is not None and file_digest(path, algorithm) != expected:
            raise ValueError(f"{algorithm.upper()} checksum does not match the built-in resource")


def valid_resource(path, resource):
    try:
        validate_resource(path, resource)
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
        if valid_resource(target, resource):
            return target
        if target.exists():
            print(f"Replacing invalid reference: {target}", file=sys.stderr)
        source = None
        candidates = bundled_candidates(kind) if builtin_paths is None else builtin_paths
        for candidate in candidates:
            candidate = Path(candidate)
            if candidate.resolve() == target.resolve():
                continue
            if valid_resource(candidate, resource):
                source = candidate
                break
            if candidate.is_file():
                print(f"Ignoring invalid bundled reference: {candidate}", file=sys.stderr)
        if source is None and resource.url is None:
            raise RuntimeError(
                f"NOT READY: the built-in {resource.classifier_id} classifier for "
                f"scikit-learn {resource.sklearn_version} has not been published, "
                f"and no compatible installed bundle was found. Install a complete "
                f"Meta2Data release containing {resource.filename}. Maintainers can "
                f"build this fixed resource with scripts/build_reference_classifiers.py "
                f"--classifier {resource.classifier_id} --output-dir <Meta2Data>/references. "
                f"The old classifier will not be downloaded or used.")
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
            validate_resource(staging, resource)
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
