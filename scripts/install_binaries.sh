#!/bin/bash
# Download prebuilt vsearch and fastp Linux binaries needed by Meta2Data.
# No compilation, no conda. Static binaries pulled from upstream releases
# and dropped into <repo>/vendor/bin alongside the pipeline source.
set -euo pipefail

VSEARCH_VERSION="2.21.1"
FASTP_VERSION="0.24.0"

# Note on vsearch version: 2.21.1 is the last release built against an older
# glibc and runs on RHEL 7 / RHEL 8 HPC nodes (glibc 2.17 / 2.28). Versions
# 2.22+ link against GLIBC_2.29 symbols and crash mid-pipeline on those
# hosts with "version `GLIBC_2.29' not found (required by vsearch)" — only
# triggered by some subcommands due to lazy symbol binding, so failures
# look intermittent. All vsearch features used by the pipeline
# (--derep_fulllength, --cluster_size, --cluster_unoise, --uchime3_denovo,
# --usearch_global, --fastq_mergepairs) are available in 2.21.1.

SCRIPT_SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_SELF_DIR")"

PREFIX="${REPO_ROOT}/vendor"
FORCE=0

usage() {
    cat <<EOF
Usage: install_binaries.sh [--prefix DIR] [--force] [-h]

Downloads vsearch ${VSEARCH_VERSION} and fastp ${FASTP_VERSION} static Linux
binaries and places them in <prefix>/bin. Idempotent unless --force is given.

Options:
  --prefix DIR   Install prefix (default: ${PREFIX})
                 The binaries land in <prefix>/bin. Meta2Data entry scripts
                 automatically prepend <repo>/vendor/bin to PATH, so leaving
                 the default means you don't have to touch your PATH.
  --force        Overwrite existing binaries even if the version matches
  -h, --help     Show this help

Supported platforms:
  vsearch: linux x86_64, linux aarch64
  fastp:   linux x86_64 only (upstream ships no aarch64 static binary)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)
            PREFIX="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown option '$1'" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "Error: Meta2Data only supports Linux; detected $(uname -s)" >&2
    exit 1
fi

ARCH="$(uname -m)"
case "$ARCH" in
    x86_64|aarch64) ;;
    *)
        echo "Error: unsupported architecture '$ARCH' (expected x86_64 or aarch64)" >&2
        exit 1
        ;;
esac

BIN_DIR="${PREFIX}/bin"
mkdir -p "$BIN_DIR"

for tool in curl tar; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "Error: '$tool' not found in PATH (required for install)" >&2
        exit 1
    }
done

TMPDIR_INSTALL="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_INSTALL"' EXIT

need_install() {
    local bin_path="$1"
    local expected_version="$2"
    [[ $FORCE -eq 1 ]] && return 0
    [[ -x "$bin_path" ]] || return 0
    local current
    current=$("$bin_path" --version 2>&1 | grep -oE "[0-9]+\.[0-9]+(\.[0-9]+)?" | head -n1 || true)
    [[ "$current" == "$expected_version" ]] && return 1
    return 0
}

install_vsearch() {
    local target="${BIN_DIR}/vsearch"
    if ! need_install "$target" "$VSEARCH_VERSION"; then
        echo "[vsearch] ${VSEARCH_VERSION} already installed at $target (use --force to reinstall)"
        return 0
    fi

    local tarball="vsearch-${VSEARCH_VERSION}-linux-${ARCH}.tar.gz"
    local url="https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VERSION}/${tarball}"
    local dest="${TMPDIR_INSTALL}/${tarball}"

    echo "[vsearch] Downloading ${url}"
    curl -fL --retry 3 --retry-delay 2 -o "$dest" "$url"

    echo "[vsearch] Extracting"
    tar -xzf "$dest" -C "$TMPDIR_INSTALL"

    local extracted="${TMPDIR_INSTALL}/vsearch-${VSEARCH_VERSION}-linux-${ARCH}/bin/vsearch"
    if [[ ! -f "$extracted" ]]; then
        echo "Error: extracted binary not found at $extracted" >&2
        return 1
    fi

    install -m 0755 "$extracted" "$target"
    echo "[vsearch] Installed to $target"
    "$target" --version 2>&1 | head -n1
}

install_fastp() {
    local target="${BIN_DIR}/fastp"
    if ! need_install "$target" "$FASTP_VERSION"; then
        echo "[fastp] ${FASTP_VERSION} already installed at $target (use --force to reinstall)"
        return 0
    fi

    if [[ "$ARCH" != "x86_64" ]]; then
        echo "Error: fastp upstream only ships a Linux x86_64 static binary." >&2
        echo "       For ${ARCH}, install fastp via your package manager or build from source." >&2
        return 1
    fi

    local url="http://opengene.org/fastp/fastp.${FASTP_VERSION}"
    local dest="${TMPDIR_INSTALL}/fastp"

    echo "[fastp] Downloading ${url}"
    curl -fL --retry 3 --retry-delay 2 -o "$dest" "$url"

    install -m 0755 "$dest" "$target"
    echo "[fastp] Installed to $target"
    "$target" --version 2>&1 | head -n1
}

echo "Installing Meta2Data runtime binaries into: $BIN_DIR"
install_vsearch
install_fastp

echo ""
echo "Done. Meta2Data entry scripts auto-prepend this directory to PATH, so"
echo "no manual PATH update is needed when you run them from '${REPO_ROOT}/bin/'."
