# Redactions

Two provenance files in this repository are redacted reading copies. Everything removed
was an absolute local filesystem path, an internal directory name, or an email address.
No hypothesis, falsification criterion, estimator specification, seed, result, or number
was changed.

## Why this file exists

The frozen analysis plan is the pre-registration. Its integrity rests on a SHA-256
recorded before any estimate was computed, and that value refers to the verbatim
original. A redacted copy cannot hash to it. Rather than publish a file whose hash
fails to verify, and rather than publish absolute paths from the author's machine, this
repository publishes the redacted copy together with this manifest, and the verbatim
original is available from the author on request.

A reader can therefore do two things. Verify this copy against the redacted hash below,
which confirms it has not been altered since release. And, on obtaining the verbatim
original, verify it against the pre-registration hash, then confirm that the only
differences are the ones enumerated here.

## Hashes

| file | SHA-256 |
|---|---|
| `results/analysis_plan.md`, verbatim original (the pre-registration) | `961cc354c62e078a49467527b35a72becabcd2020f2c03f809f69d5b3cb0b953` |
| `results/analysis_plan.md`, redacted copy as published here | `b6f7efaa0f9c7b1ee737fbb6ccec212fc327baa072c9eea6826cfc6a2a4f6c82` |
| `results/RUN_LOG.md`, verbatim original | `c467adb0e8dfb43d776109ffbc502a84e811f08a1fd1488e759a0db87a4ca5ad` |
| `results/RUN_LOG.md`, redacted copy as published here | `fbdbd6c5d6fae7ae7b23d1873a5ced550f9f19800d3c3322f90087109ab98569` |

The redacted copy of the plan also carries a provenance header that the original does
not. The header is included in its hash above.

## What was removed

### `results/analysis_plan.md`

| line | redaction |
|---|---|
| 4 | Author email address removed. The statement that the plan was prepared by an AI system on the author's behalf is retained, deliberately: it is part of the hash-frozen pre-registration and documents the preparation method from before any result existed. |
| 34 | Absolute local path to the challenge data replaced with `$ACIC22_DATA_ROOT`, the environment variable the code actually reads. |
| 37 | Absolute local path prefix removed, leaving the repository-relative path into the public benchmark repository. |
| 39 | Absolute local path to the working directory replaced with `[project working directory]`, and the reference to the source directory generalised to the data root. |

Line 26 and line 38 are unchanged. Both name repository-relative paths inside the
author's public, MIT licensed benchmark repository, which this repository cites by URL
in the README. Line 26 in particular records that the benchmark repository's Rosenbaum
module is a documented widening-bound surrogate and not a Rosenbaum randomisation
inference bound, together with the instruction that the manuscript must not describe it
as one. That is a substantive integrity constraint and is retained verbatim.

### `results/RUN_LOG.md`

| line | redaction |
|---|---|
| 122 | Internal directory names replaced with a description of what the directory holds. The scope statement itself is retained: both source trees were read only for this paper, the benchmark repository was verified unmodified, and nothing was written outside the project working directory. |
| 123 | Heading generalised and internal directory names removed. The substance is retained: only the public ACIC 2022 data files were read from the data directory, and the benchmark repository against which V1 and V2 verify is the author's own public repository. |
| 124 | One absolute local path, including a user directory and a cloud storage folder, replaced with a generic description. Two further working directory paths genericised. The technical content, namely why the packages could not be installed into the global interpreter and what was installed where, is retained in full because it is the environment disclosure. |

## Scope

Redaction was applied by exact string substitution and then verified by scanning both
output files for absolute paths, user directory names, cloud storage paths, internal
directory names, and the email address. The scan reports no remaining matches.
