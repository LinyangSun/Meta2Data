# Production and test isolation audit

Date: 19 November 2026.
Reviewed source: `a1f1ac9c191425553a119f1db7c3173ef2e4dab3`.
Comparison baseline: main `8d87dc350552c0d708e4c960bcc5bdb4537f7c92`.

## Result

No test-only computation deadline, fixture-size cap, substitute processing
script, tiny classifier, or Codespace audit path was found in the normal
Meta2Data execution path. No production-code change was required by this audit.
This finding does not mean that every timeout or numerical limit in the
repository should be removed: network safeguards, quality filters, diagnostic
sampling, and explicit validation commands have separate purposes.

The audit covered `bin/`, processing and taxonomy helpers in `scripts/`, public
parameter defaults, the environment installer and lock files, the SIF definition,
and the build workflows. It followed calls from the real entrypoints and compared
the upgrade with main, rather than treating keyword matches as defects.

## Test settings and their boundaries

| Setting | Verified location and behavior |
|---|---|
| The 30-second timeout that caused the three earlier errors | `tests/test_pip_cli.py:83–84` wraps a temporary copied test repository. That test file is unchanged from main. Normal entrypoints do not call this test wrapper. |
| Stub metadata and read-count scripts | The same unit fixture writes them below its own `TemporaryDirectory`, not over the repository's production scripts. |
| 3/5 templates and 4,000/8,000 fixture reads | Used by `tests/integration/validate_runtime.py` to generate validation inputs in its own work directory. Normal PIP/TAXA does not invoke that generator. |
| Two threads and one concurrent dataset used in the Codespace audit | Arguments to individual audit invocations. Public defaults remain four total threads and two concurrent PIP datasets; users can override them. |
| Confidence 0.8 used in the cache regression | An explicit argument in the audit driver. The public TAXA default is still 0.7 and normal classification receives the requested value. |
| Historical Codespace paths and audit drivers | Present in documentation/evidence only. They are not in the normal runtime PATH and are not imported or sourced by production entrypoints. |
| Tiny classifier used for an earlier training roundtrip | Not registered or selected by normal resource preparation. Production validates the chosen model's semantic type and scikit-learn version, plus the published SILVA 144 identity/checksum. |
| SILVA 144 rank-schema assertions | Used by the integration validator. Production does not truncate taxonomy strings or require every call to reach genus. |

`bin/Meta2Data-AmpliconPIP`, `scripts/run.sh`, `scripts/parameters.py`, and
`docs/parameters.default.json` are byte-identical to main. The PIP functions file
differs only by four required DADA2 base-transition output arguments.
TAXA supplies the complete oriented sequence artifact to classification and
passes the user's CPU and confidence settings. Input collection has no test
dataset-count cap.

## Existing features retained intentionally

- **Explicit `--test`:** this pre-existing public option remains off by default.
  Only opting into it activates the two-rows-per-project metadata subset. Normal
  runs do not enter that branch.
- **Diagnostic sampling:** primer/quality checks and length/N-statistic estimates
  inspect subsets to choose processing parameters. The subsequent processing
  loops still consume the full applicable input files. These existing sampling
  settings are not the validation fixture's read-count limits.
- **Network timeouts:** HTTP requests and reference connection setup retain
  bounded waits. They do not impose a total QIIME or FASTQ-processing deadline.
- **Installation probes:** the SRA installer allows 30 seconds for `--version`
  and `--help` checks. This timeout is not applied to SRA extraction or analysis.
- **Existing tool/resource policies:** `fasterq-dump --threads 2` predates this
  upgrade. The maintainer-only SILVA 138.1 preparation recipe also has explicit
  single-thread preparation settings; normal analysis never invokes it.
- **Metadata validation:** the taxonomy cache checker caps each small YAML
  metadata/action record at 1 MiB. It does not cap model/archive size, reads,
  features, samples, or taxonomy depth.

## Native and SIF execution

Native installation adds the repository's `bin/` directory to PATH. The SIF
runscript directly executes the real Meta2Data entrypoint with the user's
arguments. Neither inserts fixture limits, test wrappers, or audit paths.

The SIF bundles test and documentation files for validation and reference, but
ordinary analysis does not execute them. `%test`, `--smoke`, and the GitHub
Actions time limits apply to build/verification jobs only. The release validator
already runs unit tests before real workflows, sequentially; the earlier ad hoc
Codespace overlap is not part of that CI order.

## Direct execution check

In an isolated local audit directory, the unchanged production resource profiler
ran a child command that waited 31 seconds and then printed a completion marker.
The wrapper returned **0** after **31.1149 seconds**; its own event recorded
**31.0247 seconds** and success. This directly checks that the profiler does not
apply the test suite's 30-second cutoff. It is a local wrapper probe, not an
additional QIIME workflow run.

The profiler waits for the child with `os.wait4(child.pid, 0)`. Its separate
five-second monitoring-thread join happens after that child has exited.
The exact command and result are preserved in
[production-timeout-probe.json](production-timeout-probe.json).

No classifier training, package rebuilding, or main-branch changes were made.
