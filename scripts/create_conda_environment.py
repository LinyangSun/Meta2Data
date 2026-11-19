#!/usr/bin/env python3
"""Install Meta2Data's locked Linux runtime and pinned standalone SRA Toolkit.

env.yml lists direct requirements. Release constraints retain compatible
dependency versions/builds; the explicit lock is tied to both files by SHA256.
Configured extra channels and default packages cannot change this installation.
No Conda configuration is changed.
"""
import argparse
import hashlib
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
PIN = re.compile(r"([A-Za-z0-9_.-]+)=[A-Za-z0-9_.+!-]+\Z")


def read_freeze(path):
    """Read this project's deliberately flat YAML subset, rejecting other forms."""
    sections = {}
    section = None
    for number, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("name:") and section is None:
            continue
        if line in ("channels:", "dependencies:"):
            section = line[:-1]
            if section in sections:
                raise ValueError(f"Duplicate {section} section in {path}")
            sections[section] = []
        elif line.startswith("- ") and section is not None:
            sections[section].append(line[2:].strip())
        else:
            raise ValueError(f"Unsupported env.yml structure at line {number}: {line}")
    if set(sections) != {"channels", "dependencies"}:
        raise ValueError("env.yml must contain channels and dependencies")
    channels = [item for item in sections["channels"] if item != "nodefaults"]
    if not channels or "defaults" in channels or len(channels) != len(set(channels)):
        raise ValueError("env.yml must specify unique non-default channels")
    pins = sections["dependencies"]
    if not pins or any(PIN.fullmatch(pin) is None for pin in pins):
        raise ValueError("env.yml dependencies must be flat, single-version Conda pins")
    names = [pin.split("=", 1)[0] for pin in pins]
    if len(names) != len(set(names)):
        raise ValueError("env.yml contains duplicate package requirements")
    return channels, pins


def read_constraints(path):
    packages = {}
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("=")
        if len(parts) != 3 or not all(re.fullmatch(r"[A-Za-z0-9_.+!-]+", part)
                                      for part in parts):
            raise ValueError(f"Release constraint must specify name=version=build: {line}")
        name, version, build = parts
        if name in packages:
            raise ValueError(f"Duplicate release constraint: {name}")
        packages[name] = (version, build)
    if not packages:
        raise ValueError("Release constraints must not be empty")
    return packages


def read_lock(path, freeze, channels, pins, constraints=None):
    """Reject stale, incomplete, or foreign-channel explicit package lists."""
    lines = Path(path).read_text().splitlines()
    digest = hashlib.sha256(Path(freeze).read_bytes()).hexdigest()
    if f"# env.yml sha256: {digest}" not in lines:
        raise ValueError("The environment lock does not match env.yml; regenerate both together")
    constraints = Path(constraints or ROOT / "environment-constraints-linux-64.txt")
    digest = hashlib.sha256(constraints.read_bytes()).hexdigest()
    if f"# release constraints sha256: {digest}" not in lines:
        raise ValueError("The environment lock does not match release constraints")
    expected = read_constraints(constraints)
    specifications = [line for line in lines if line and not line.startswith("#")]
    if not specifications or specifications.pop(0) != "@EXPLICIT":
        raise ValueError("The environment lock must be an explicit package list")
    bases = [channel if channel.startswith("https://") else
             f"https://conda.anaconda.org/{channel}" for channel in channels]
    packages = {}
    for line in specifications:
        url, separator, digest = line.rpartition("#")
        if not separator or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError("Every locked package must include its SHA256")
        source = urlsplit(url)
        filename = source.path.rsplit("/", 1)[-1]
        allowed = [f"{base}/{subdir}/{filename}" for base in bases
                   for subdir in ("linux-64", "noarch")]
        if url not in allowed:
            raise ValueError(f"Unsupported package source in lock: {url}")
        suffix = next((suffix for suffix in (".conda", ".tar.bz2")
                       if filename.endswith(suffix)), None)
        if suffix is None:
            raise ValueError(f"Unsupported package archive: {filename}")
        try:
            name, version, build = filename[:-len(suffix)].rsplit("-", 2)
        except ValueError as error:
            raise ValueError(f"Malformed package filename: {filename}") from error
        if name in packages:
            raise ValueError(f"Duplicate package in lock: {name}")
        packages[name] = (version, build)
    for pin in pins:
        name, version = pin.split("=")
        if packages.get(name, (None, None))[0] != version:
            raise ValueError(f"The environment lock does not provide {pin}")
    if packages != expected:
        raise ValueError("The environment lock differs from the exact release package constraints")
    return specifications


def create(prefix, conda, *, dry_run=False, freeze=None, lock=None):
    freeze, lock = Path(freeze or ROOT / "env.yml"), Path(lock or ROOT / "environment-linux-64.lock")
    channels, pins = read_freeze(freeze)
    read_lock(lock, freeze, channels, pins)
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "amd64"):
        raise ValueError("This runtime supports Linux x86_64 only")
    libc = os.confstr("CS_GNU_LIBC_VERSION") or ""
    match = re.fullmatch(r"glibc (\d+)\.(\d+)", libc)
    if match is None or tuple(map(int, match.groups())) < (2, 28):
        raise ValueError(f"This runtime requires glibc >=2.28; detected {libc or 'unknown libc'}")
    prefix = Path(prefix).expanduser().absolute()
    # conda create --yes can replace an existing environment. Never allow that.
    if prefix.is_symlink() or (prefix.exists() and
                              (not prefix.is_dir() or any(prefix.iterdir()))):
        raise ValueError(f"Destination must be absent or empty: {prefix}")
    executable = shutil.which(conda)
    if executable is None:
        raise ValueError(f"Conda executable not found: {conda}; supply --conda /path/to/conda")
    command = [executable, "create", "--yes", "--solver", "libmamba",
               "--prefix", str(prefix), "--override-channels", "--no-default-packages"]
    for channel in channels:
        command.extend(["--channel", channel])
    command.extend(["--file", str(lock)])
    if dry_run:
        command.append("--dry-run")
    print(shlex.join(command), flush=True)
    result = subprocess.run(
        command, env=dict(os.environ, CONDA_CHANNEL_PRIORITY="flexible"))
    if result.returncode or dry_run:
        return result.returncode
    return subprocess.run([str(prefix / "bin/python"), str(ROOT / "scripts/install_sra_toolkit.py"),
                           "--prefix", str(prefix)]).returncode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True,
                        help="New or empty destination; existing environments are never replaced")
    parser.add_argument("--conda", default=os.environ.get("CONDA_EXE") or "conda",
                        help="Conda executable (default: CONDA_EXE or conda on PATH)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Check the lock and Conda command without downloading/installing packages")
    args = parser.parse_args(argv)
    try:
        return create(args.prefix, args.conda, dry_run=args.dry_run)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Environment creation failed: {error}\n")


if __name__ == "__main__":
    sys.exit(main())
