#!/bin/bash
# Pre-flight dependency check for Meta2Data.
# Can be run directly, or sourced by an entry script which then calls
# meta2data_check_dependencies or meta2data_ensure_python_deps.
#
# Exits (or returns when sourced) non-zero if any required binary, QIIME2
# plugin, or Python package is missing. Prints a report with install hints.
#
# meta2data_ensure_python_deps is a lighter-weight helper used by each
# entry script at startup: it checks the Python packages Meta2Data needs,
# and auto-installs any missing ones into the active python3 environment
# via `python3 -m pip install`. Set META2DATA_SKIP_DEP_CHECK=1 to bypass.
#
# meta2data_ensure_vendor_binaries is the parallel helper for vsearch and
# fastp. AmpliconPIP / AmpliconTAXA entry scripts call it at startup; it adds
# <repo>/vendor/bin to PATH and, if either binary is still missing, fails with
# install hints (conda / module). MetaDL does not need these binaries and does
# not call this helper. Set META2DATA_SKIP_DEP_CHECK=1 to bypass.

# Required non-QIIME2 binaries (vsearch / fastp).
_M2D_VENDOR_BINARIES=(
    vsearch
    fastp
)

# Long-read (ONT) tools — needed ONLY by the OXFORD_NANOPORE branch in run.sh.
# Reported as an advisory in the dependency check, but NEVER used to gate other
# platforms: run.sh checks these lazily, per ONT dataset, so an Illumina/454/
# Ion Torrent/PacBio run is never blocked by a missing long-read tool.
_M2D_ONT_BINARIES=(
    chopper
    minimap2
    racon
)

# Required generic system utilities (expected on any Linux distro).
_M2D_SYSTEM_BINARIES=(
    python3
    wget
    curl
    unzip
    tar
    gzip
    md5sum
    awk
    sed
    grep
    find
)

# Required QIIME2 plugins. Names match the strings `qiime info` prints.
_M2D_QIIME_PLUGINS=(
    cutadapt
    dada2
    demux
    feature-classifier
    feature-table
    fragment-insertion
    quality-filter
    rescript
    vsearch
)

# Python packages NOT auto-provided by QIIME2 (biopython) plus a sanity
# check on ones that ARE (pandas/numpy/requests/biom-format) to catch a
# broken QIIME2 install.
#
# Each entry is "import_name:pip_name".
_M2D_PYTHON_PACKAGES=(
    "Bio:biopython"
    "pandas:pandas"
    "numpy:numpy"
    "requests:requests"
    "biom:biom-format"
)

_m2d_missing_vendor=()
_m2d_missing_system=()
_m2d_missing_qiime=()
_m2d_missing_python=()

_m2d_check_binary() {
    command -v "$1" >/dev/null 2>&1
}

_m2d_check_python_package() {
    python3 -c "import $1" >/dev/null 2>&1
}

# Collect QIIME2 plugin list once and grep locally — avoids launching
# `qiime <plugin> --help` N times (each spawn costs ~1s).
_m2d_load_qiime_plugin_list() {
    if ! command -v qiime >/dev/null 2>&1; then
        _M2D_QIIME_INFO=""
        return 1
    fi
    _M2D_QIIME_INFO=$(qiime info 2>&1 || true)
    return 0
}

_m2d_check_qiime_plugin() {
    grep -qE "^[[:space:]]*$1:[[:space:]]" <<<"$_M2D_QIIME_INFO"
}

meta2data_check_dependencies() {
    _m2d_missing_vendor=()
    _m2d_missing_system=()
    _m2d_missing_qiime=()
    _m2d_missing_python=()

    for bin in "${_M2D_VENDOR_BINARIES[@]}"; do
        _m2d_check_binary "$bin" || _m2d_missing_vendor+=("$bin")
    done

    for bin in "${_M2D_SYSTEM_BINARIES[@]}"; do
        _m2d_check_binary "$bin" || _m2d_missing_system+=("$bin")
    done

    if command -v qiime >/dev/null 2>&1; then
        _m2d_load_qiime_plugin_list || true
        for plugin in "${_M2D_QIIME_PLUGINS[@]}"; do
            _m2d_check_qiime_plugin "$plugin" || _m2d_missing_qiime+=("$plugin")
        done
    else
        _m2d_missing_qiime=("${_M2D_QIIME_PLUGINS[@]}")
        _m2d_missing_system+=("qiime")
    fi

    for entry in "${_M2D_PYTHON_PACKAGES[@]}"; do
        local import_name="${entry%%:*}"
        local pip_name="${entry##*:}"
        if command -v python3 >/dev/null 2>&1; then
            _m2d_check_python_package "$import_name" || _m2d_missing_python+=("$pip_name")
        else
            _m2d_missing_python+=("$pip_name")
        fi
    done

    # Advisory only: long-read (ONT) tools are needed solely by the
    # OXFORD_NANOPORE branch, so they are reported but never counted toward the
    # failure total. run.sh checks them lazily per ONT dataset.
    local _m2d_missing_ont=()
    for bin in "${_M2D_ONT_BINARIES[@]}"; do
        _m2d_check_binary "$bin" || _m2d_missing_ont+=("$bin")
    done
    if [[ ${#_m2d_missing_ont[@]} -gt 0 ]]; then
        echo "[check] Note: ONT-only tools not found: ${_m2d_missing_ont[*]}"
        echo "        (needed only for Oxford Nanopore datasets; install:"
        echo "         conda install -c bioconda ${_m2d_missing_ont[*]})"
    fi

    local total=$(( ${#_m2d_missing_vendor[@]} + ${#_m2d_missing_system[@]} + ${#_m2d_missing_qiime[@]} + ${#_m2d_missing_python[@]} ))

    if [[ $total -eq 0 ]]; then
        echo "[check] All dependencies satisfied."
        return 0
    fi

    echo ""
    echo "=============================================================="
    echo "  Meta2Data dependency check FAILED — ${total} item(s) missing"
    echo "=============================================================="

    if [[ ${#_m2d_missing_vendor[@]} -gt 0 ]]; then
        echo ""
        echo "Missing pipeline binaries:"
        for b in "${_m2d_missing_vendor[@]}"; do
            echo "  - $b"
        done
        echo "  Install: conda install -c bioconda ${_m2d_missing_vendor[*]}"
    fi

    if [[ ${#_m2d_missing_system[@]} -gt 0 ]]; then
        echo ""
        echo "Missing system utilities (should be available on any Linux distro):"
        for b in "${_m2d_missing_system[@]}"; do
            echo "  - $b"
        done
        echo "  Install via your package manager, e.g.:"
        echo "    sudo apt-get install ${_m2d_missing_system[*]}"
    fi

    if [[ ${#_m2d_missing_qiime[@]} -gt 0 ]]; then
        echo ""
        if [[ -z "${_M2D_QIIME_INFO:-}" ]]; then
            echo "QIIME2 not detected at all. Install the 2024.10 amplicon distribution"
            echo "following https://docs.qiime2.org/2024.10/install/, then activate its"
            echo "conda environment before running Meta2Data."
        else
            echo "Missing QIIME2 plugins:"
            for p in "${_m2d_missing_qiime[@]}"; do
                echo "  - $p"
            done
            echo "  Your QIIME2 env looks incomplete. Install the full amplicon"
            echo "  distribution (not just the core) per QIIME2's official docs."
        fi
    fi

    if [[ ${#_m2d_missing_python[@]} -gt 0 ]]; then
        echo ""
        echo "Missing Python packages:"
        for p in "${_m2d_missing_python[@]}"; do
            echo "  - $p"
        done
        echo "  Install with:"
        echo "    pip install ${_m2d_missing_python[*]}"
    fi

    echo ""
    echo "=============================================================="
    return 1
}

# -----------------------------------------------------------------------------
# Lightweight Python-dep ensure helper (called by every Meta2Data entry script)
# -----------------------------------------------------------------------------
# Auto-installs any missing Python packages into whatever python3 is currently
# on PATH. Uses `python3 -m pip install`, so the install lands in the exact
# site-packages the import check runs against (system python / venv / conda
# env — Meta2Data does not care which, it just follows the active interpreter).
#
# Design choices:
#   - Unified dep list across all subcommands. MetaDL, AmpliconPIP and AmpliconTAXA
#     all end up touching one of biopython / pandas / numpy / requests, so
#     there is no point in per-subcommand grouping.
#   - No interactive prompt. If the user ran a Meta2Data command, they've
#     opted in to whatever it needs; blocking on stdin breaks CI and HPC jobs.
#   - Graceful bypass via META2DATA_SKIP_DEP_CHECK=1 for users who manage their
#     environment by hand and want Meta2Data to stay out of pip.

_M2D_REQUIRED_PYDEPS=(
    "Bio:biopython"
    "pandas:pandas"
    "numpy:numpy"
    "requests:requests"
)

meta2data_ensure_python_deps() {
    [[ "${META2DATA_SKIP_DEP_CHECK:-}" == "1" ]] && return 0

    if ! command -v python3 >/dev/null 2>&1; then
        echo "Error: python3 not found in PATH." >&2
        echo "       Meta2Data needs a Python 3 interpreter on PATH." >&2
        echo "       Install python3 via your package manager or activate a" >&2
        echo "       conda/venv environment that provides it, then retry." >&2
        return 1
    fi

    local -a missing=()
    local entry import_name pip_name
    for entry in "${_M2D_REQUIRED_PYDEPS[@]}"; do
        import_name="${entry%%:*}"
        pip_name="${entry##*:}"
        python3 -c "import ${import_name}" >/dev/null 2>&1 \
            || missing+=("$pip_name")
    done

    [[ ${#missing[@]} -eq 0 ]] && return 0

    local py_prefix
    py_prefix=$(python3 -c 'import sys; print(sys.prefix)' 2>/dev/null || echo "unknown")
    echo "[deps] Missing Python packages: ${missing[*]}"
    echo "[deps] Installing into active python3 env: ${py_prefix}"
    if ! python3 -m pip install "${missing[@]}"; then
        echo "" >&2
        echo "Error: failed to auto-install Python packages: ${missing[*]}" >&2
        echo "       This can happen on read-only environments or when pip is" >&2
        echo "       unavailable. Install them manually with:" >&2
        echo "         python3 -m pip install ${missing[*]}" >&2
        echo "       Or set META2DATA_SKIP_DEP_CHECK=1 to bypass this check." >&2
        return 1
    fi
    echo "[deps] OK"
    return 0
}

# -----------------------------------------------------------------------------
# Lightweight vendor-binary ensure helper
# -----------------------------------------------------------------------------
# Ensures vsearch / fastp are on PATH for the pipeline. Prepends
# <repo>/vendor/bin to PATH and, if either binary is still missing, fails with
# install hints (conda / module). Called at startup by AmpliconPIP /
# AmpliconTAXA entry scripts only — MetaDL does not need these.

meta2data_ensure_vendor_binaries() {
    [[ "${META2DATA_SKIP_DEP_CHECK:-}" == "1" ]] && return 0

    local scripts_dir repo_root vendor_bin
    scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    repo_root="$(dirname "$scripts_dir")"
    vendor_bin="${repo_root}/vendor/bin"

    # Prepend vendor/bin to PATH so command -v can find binaries there too.
    [[ -d "$vendor_bin" ]] && export PATH="${vendor_bin}:${PATH}"

    local missing=()
    local b
    for b in "${_M2D_VENDOR_BINARIES[@]}"; do
        command -v "$b" >/dev/null 2>&1 || missing+=("$b")
    done

    [[ ${#missing[@]} -eq 0 ]] && return 0

    echo "" >&2
    echo "Error: missing required binaries: ${missing[*]}" >&2
    echo "       Install them via one of:" >&2
    echo "         - conda install -c bioconda ${missing[*]}" >&2
    echo "         - module load ${missing[*]}" >&2
    echo "       Or set META2DATA_SKIP_DEP_CHECK=1 to bypass this check." >&2
    return 1
}

# When executed directly, run the check and exit with its status.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    meta2data_check_dependencies
    exit $?
fi
