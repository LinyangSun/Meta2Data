# Codespace compatibility checks, 19 November 2026

Tested source: b6d7bd8d7241289ed618917887349867cbadfdad.
Environment: QIIME 2 2026.7.0 and scikit-learn 1.7.1, installed from the
485-package minimal release lock plus the pinned standalone SRA Toolkit.

## Observed results

- Runtime interfaces, executable checks, entrypoint help, and all 304 unit tests passed. No tests were skipped.
- The original main-branch GG2 classifier was tested directly through qiime feature-classifier classify-sklearn, without changing its metadata or the QIIME compatibility check. It exited 1 before classification: the artifact records scikit-learn 1.4.2, while the environment uses 1.7.1. No taxonomy artifact was produced.
- The official full-length SILVA 144 classifier annotated the three supplied representative sequences. Its observed prefixes were d__;k__;p__;c__;o__;f__;g__. This confirms the tested output includes Kingdom between Domain and Phylum.
- A separate two-sequence synthetic classifier passed training, saving, reloading, and classification in the new environment.

## Evidence

- runtime-checks.json: actual environment and unit-test results.
- legacy-gg2-compatibility.json: original GG2 URL, UUID, SHA256, exact command, and exit status.
- legacy-gg2-classify.log: unedited QIIME exception and traceback.
- classifier-probes.json: tiny training roundtrip and official SILVA 144 model identity/output.
- silva144-taxonomy.tsv and representative-sequences.fasta: the small classification input and output.

These checks do not establish biological accuracy or complete release validation.
The full platform/workflow matrix, compatible replacement legacy models, and
fresh SIF validation are still pending. No full reference classifier was trained.
