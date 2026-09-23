# Profiling is opt-in. PATH wrappers preserve stdin/stdout, argv and exit status,
# and also observe supported tools launched by Python helpers.
# M2D_PROFILE_STAGE labels one captured invocation only; SAMPLE and
# DECISION_REPORT persist to its children. Invoke Python helpers through PATH
# (python3), so their own stage and child hierarchy are both captured.
if [[ -n "${M2D_PROFILE_DIR:-}" && "${M2D_PROFILE_ACTIVE:-0}" != 1 ]]; then
    [[ "$M2D_PROFILE_DIR" = /* ]] || { echo 'M2D_PROFILE_DIR must be absolute' >&2; exit 2; }
    export M2D_PROFILE_PYTHON="$(python3 -c 'import sys; print(sys.executable)')"
    _profile_bin="${M2D_PROFILE_DIR}/.bin-${BASHPID}"
    python3 "${SCRIPTS}/resource_profile.py" bootstrap --directory "$_profile_bin" || exit 2
    export PATH="${_profile_bin}:${PATH}" M2D_PROFILE_ACTIVE=1
    unset _profile_bin
fi
