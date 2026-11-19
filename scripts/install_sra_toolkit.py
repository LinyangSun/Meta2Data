#!/usr/bin/env python3
"""Install the pinned upstream SRA Toolkit beside a Linux x86_64 Conda runtime.

The official binary bundle avoids the conflicting libuuid/ossuuid files in the
Bioconda package. Keep its dispatcher, original executables, and schemas
together; expose only the two commands used by Meta2Data.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import posixpath
import re
import subprocess
import tarfile
import tempfile
import urllib.request


VERSION = "3.4.1"
DIRECTORY = f"sratoolkit.{VERSION}-alma_linux64"
ARCHIVE_URL = f"https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/{VERSION}/{DIRECTORY}.tar.gz"
ARCHIVE_SHA256 = "f4850a32975596ca0755d928b166ecb0095918284ea992fba10f1acb3fcff64d"
ARCHIVE_MD5 = "3758b10c272969852e6fd827739eb93b"
ARCHIVE_BYTES = 89143374
MANIFEST = "meta2data-install.json"
# Verified from the checksum-pinned upstream archive. These three files form
# the dispatcher and the actual commands needed by the pipeline.
PAYLOAD_SHA256 = {
    "bin/sratools.3.4.1": "4dd8a1e0996aef1d4879632da88908351503528d653d6fc6e1bc98ceabc8bd40",
    "bin/prefetch-orig.3.4.1": "e4d7479c1e912005e82818decbe47181be0a35dce9fdd3284ba6516ee49abebb",
    "bin/fasterq-dump-orig.3.4.1": "b3a4dcd4078c1867bc980fecde4333efadd435e5867d7a834cd3b0a40c419d9f",
}
TOOL_FLAGS = {
    "prefetch": ("--max-size", "--output-directory"),
    "fasterq-dump": ("--split-files", "--skip-technical", "--threads", "--outdir"),
}


def digest(path, algorithm="sha256"):
    result = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def check_platform():
    if platform.system() != "Linux" or platform.machine() != "x86_64":
        raise RuntimeError("This SRA Toolkit bundle supports Linux x86_64 only")
    libc, version = platform.libc_ver()
    if libc != "glibc" or tuple(map(int, version.split(".")[:2])) < (2, 28):
        raise RuntimeError("The pinned SRA Toolkit requires glibc 2.28 or later")


def download(destination):
    with urllib.request.urlopen(ARCHIVE_URL, timeout=60) as response, \
            Path(destination).open("wb") as stream:
        for block in iter(lambda: response.read(1024 * 1024), b""):
            stream.write(block)


def validate_download(path):
    if Path(path).stat().st_size != ARCHIVE_BYTES:
        raise ValueError("SRA Toolkit archive has an unexpected size")
    if digest(path) != ARCHIVE_SHA256 or digest(path, "md5") != ARCHIVE_MD5:
        raise ValueError("SRA Toolkit archive checksum does not match the pinned release")


def extract_archive(path, destination):
    """Accept only this bundle's root and links contained inside that root."""
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        seen = set()
        for member in members:
            name = PurePosixPath(member.name)
            if (name.is_absolute() or not name.parts or name.parts[0] != DIRECTORY
                    or ".." in name.parts or "\\" in member.name):
                raise ValueError(f"Unsafe SRA archive path: {member.name}")
            normalized = name.as_posix()
            if normalized in seen:
                raise ValueError(f"Duplicate SRA archive path: {member.name}")
            seen.add(normalized)
            if not (member.isfile() or member.isdir() or member.issym()):
                raise ValueError(f"Unsupported SRA archive entry: {member.name}")
            if member.issym():
                target = PurePosixPath(member.linkname)
                resolved = posixpath.normpath(str(name.parent / target))
                if (target.is_absolute() or "\\" in member.linkname
                        or not resolved.startswith(DIRECTORY + "/")):
                    raise ValueError(f"SRA archive link escapes its directory: {member.name}")
        if not hasattr(tarfile, "data_filter"):
            raise RuntimeError("A Python runtime with safe tarfile extraction is required")
        archive.extractall(destination, members=members, filter="data")


def payload_hashes(root):
    values = {}
    for path in sorted(root.rglob("*")):
        if path.is_file() and not path.is_symlink() and path.name != MANIFEST:
            values[path.relative_to(root).as_posix()] = digest(path)
    for name, expected in PAYLOAD_SHA256.items():
        if values.get(name) != expected:
            raise ValueError(f"SRA Toolkit executable checksum mismatch: {name}")
    if not any(name.startswith("schema/") for name in values):
        raise ValueError("SRA Toolkit schema files are missing")
    return values


def payload_symlinks(root):
    values = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            if not path.resolve(strict=True).is_relative_to(root.resolve()):
                raise ValueError(f"SRA Toolkit link escapes its installation: {path}")
            values[path.relative_to(root).as_posix()] = os.readlink(path)
    return values


def verify_commands(bin_directory, runner=None):
    runner = runner or subprocess.run
    versions = {}
    for tool, flags in TOOL_FLAGS.items():
        executable = bin_directory / tool
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError(f"SRA Toolkit command is missing or not executable: {executable}")
        for option in ("--version", "--help"):
            result = runner([str(executable), option], check=True, text=True,
                            capture_output=True, timeout=30)
            output = result.stdout + result.stderr
            if option == "--version":
                match = re.search(r":\s*(\d+\.\d+\.\d+)\s*$", output)
                if not match or match.group(1) != VERSION:
                    raise ValueError(f"Unexpected SRA Toolkit version: {output.strip()}")
                versions[tool] = match.group(1)
            elif any(flag not in output for flag in flags):
                raise ValueError(f"{tool} is missing a command-line option required by Meta2Data")
    return versions


def install_toolkit(prefix, *, downloader=None, runner=None):
    check_platform()
    prefix = Path(prefix).resolve()
    if not prefix.is_dir():
        raise ValueError("--prefix must be an existing Conda environment directory")
    opt, binary = prefix / "opt", prefix / "bin"
    opt.mkdir(exist_ok=True)
    binary.mkdir(exist_ok=True)
    target = opt / DIRECTORY
    manifest_path = target / MANIFEST
    links = {name: Path("../opt") / DIRECTORY / "bin" / name for name in TOOL_FLAGS}
    with (opt / ".meta2data-sra.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for name, relative in links.items():
            entry = binary / name
            if os.path.lexists(entry) and not (
                    entry.is_symlink() and os.readlink(entry) == str(relative)):
                raise FileExistsError(f"Refusing to overwrite existing command: {entry}")
        if os.path.lexists(target):
            if target.is_symlink() or not manifest_path.is_file():
                raise FileExistsError(f"Refusing to replace an unregistered toolkit: {target}")
            manifest = json.loads(manifest_path.read_text())
            if (manifest.get("version") != VERSION or
                    manifest.get("archive", {}).get("sha256") != ARCHIVE_SHA256):
                raise ValueError(f"Installed SRA Toolkit provenance differs from this release: {target}")
            if payload_hashes(target) != manifest.get("payload_sha256"):
                raise ValueError(f"Installed SRA Toolkit files changed: {target}")
            if payload_symlinks(target) != manifest.get("payload_symlinks"):
                raise ValueError(f"Installed SRA Toolkit links changed: {target}")
        else:
            with tempfile.TemporaryDirectory(prefix=".sra-toolkit-", dir=opt) as temporary:
                staging = Path(temporary)
                archive = staging / "toolkit.tar.gz"
                (downloader or download)(archive)
                validate_download(archive)
                extract_archive(archive, staging)
                extracted = staging / DIRECTORY
                hashes = payload_hashes(extracted)
                symlinks = payload_symlinks(extracted)
                versions = verify_commands(extracted / "bin", runner)
                manifest = {
                    "schema_version": 1, "version": VERSION, "platform": "linux-x86_64",
                    "installed_utc": datetime.now(timezone.utc).isoformat(),
                    "archive": {"url": ARCHIVE_URL, "sha256": ARCHIVE_SHA256,
                                "md5": ARCHIVE_MD5, "bytes": ARCHIVE_BYTES},
                    "tools": {name: {"version": version, "path": f"bin/{name}"}
                              for name, version in versions.items()},
                    "payload_sha256": hashes,
                    "payload_symlinks": symlinks,
                }
                (extracted / MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
                os.rename(extracted, target)
        for name, relative in links.items():
            entry = binary / name
            if not os.path.lexists(entry):
                entry.symlink_to(relative)
        # Exercise the actual exposed paths to catch dispatcher lookup errors.
        verify_commands(binary, runner)
    return manifest_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = install_toolkit(args.prefix)
    except (OSError, ValueError, RuntimeError, tarfile.TarError,
            subprocess.SubprocessError) as error:
        parser.exit(1, f"SRA Toolkit installation failed: {error}\n")
    print(f"SRA Toolkit {VERSION} ready; provenance: {manifest}")


if __name__ == "__main__":
    main()
