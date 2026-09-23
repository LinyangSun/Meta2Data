# Revision validation

The revision was implemented in `Meta2Data-main2`. The source tree before edits
was copied to `/tmp/meta2data-main2-before-revision` for comparison and rollback.

Run the portable regression suite from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

Final result after the code and language review: **88 tests passed**.

The suite covers common paired filename forms, lanes/chunks, source preservation,
exact sample IDs, empty/missing input diagnostics, b1 decision boundaries,
unknown-primer skipping, mixed orientations, configuration precedence, restart
invalidation, recursive artifact discovery/deduplication, and all four
method/region TAXA combinations, and sequence-derived vsearch feature IDs with
unchanged sample count totals. The follow-up review also covers source paths and
symlinks overlapping managed output, local `-m .`, actual R2 decisions in mixed
orientations, explicit-primer tool failures, legitimate `final-*` dataset names,
required TAXA failure propagation, and independent GG2/SILVA caches. Public-entrypoint tests use isolated package
copies and substitute the computational tools; they do not claim biological
pipeline validation.

An independent comparison checked 860 deterministic detector cases against
`../Meta2Data-main/scripts/entropy_primer_detect.py`, changing only its fold
threshold from 32 to the b1 value 16. The cases comprise 552 known-oligo cases
(all 23 bundled primers, IUPAC expansions, orientations, offsets and richness),
220 unprimed mixtures, 48 fold boundaries and 40 supported-frequency boundaries.
There were zero behavioral differences in detection or trim length. For the
327 database matches, an uncomputed fold is now recorded as JSON `null` rather
than the old placeholder zero. This is a detector parity audit, not a rerun of
the original b1 simulation. The audit was repeated after the mixed-PE correction
and retained zero behavioral differences across all 860 cases.

All seven shell files were checked with `bash -n`, and all 20 Python sources
were parsed. All four public entry points returned successful English `--help`
output. A scan of 40 source, configuration and documentation files found no
literal Chinese, Japanese or Korean text. Localized field labels needed to
parse external CNCB HTML remain supported through Unicode escapes in regexes;
logs and documentation are in English. QIIME2, fastp, vsearch and the scientific Python dependencies were not
available in the active shell, so no real denoising, database download or complete
biological end-to-end run was performed in this session.
