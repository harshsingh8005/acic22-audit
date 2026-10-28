# A Pre-Registered Ground-Truth Audit of Interval Coverage and Sensitivity Analysis on ACIC 2022

Analysis code for the pre-registered audit of confidence interval coverage and of a
unified sensitivity analysis layer (Rosenbaum bounds, E-values, omitted variable bias
bounds) across all 3,400 practice level cohorts of the ACIC 2022 Track 2 data
challenge, where the true sample average treatment effects on the treated are public.

Nothing in this repository is a sensitivity analysis method proposal. The code
evaluates existing methods against a known truth under a controlled omission.

## Pre-registration

The analysis plan was frozen and hash logged before any estimate was computed.

| | |
|---|---|
| Frozen plan | `results/analysis_plan.md` |
| SHA-256 | `961cc354c62e078a49467527b35a72becabcd2020f2c03f809f69d5b3cb0b953` |
| Frozen at | 2026-08-18T20:52:50Z |

`results/RUN_LOG.md` is the execution record: every command, wall time, output file
with its SHA-256, and every deviation from the plan, in the order they happened.
Read it alongside this README. Where the two disagree, the run log is authoritative.

Both files are published here as redacted reading copies. The SHA-256 above refers to
the verbatim original of the plan, which is available from the author on request; the
copy in this repository has its own hash and cannot match that value. `REDACTIONS.md`
lists every removal line by line, gives both hashes for both files, and records that
only absolute local filesystem paths, internal directory names and one email address
were removed. No hypothesis, criterion, specification, seed or number differs.

## Layout

```
code/                    17 analysis files, ~2,850 lines
  acic_common.py         cohort loading, DiD reduction, folds, stratified sampling
  estimators.py          E0 OLS-ATT, E1 AIPW-ATT, E2 DR-ATT with Crump trim, E6 ridge T-learner
  sensitivity.py         Rosenbaum bounds, E-values, OVB robustness values and bias bounds
  neural.py              TARNet and DragonNet (PyTorch)
  validate_software.py   software checks V1 to V6
  run_fast_arm.py        confirmatory arm, 3,400 cohorts
  run_neural_arm.py      neural arm, 170 cohorts
  run_ovb_grid.py        OVB coverage grid for figure F4
  run_e3_e8_fast.py      E3 DoubleML, E8 EconML            (cross library environment)
  run_catenets.py        catenets replication of E4 and E5  (cross library environment)
  run_e7_dowhy.py        E7 DoWhy cross check
  analyze.py             hypothesis tests, strata summaries, calibration curve
  analyze_addendum.py    cross library analysis
  figures.py             F1 to F6
  figure_fidelity_check.py, assumption_checks_report.py, make_report_tables.py
results/
  analysis_plan.md       the frozen pre-registration
  RUN_LOG.md             execution record with output hashes
  tables/, figures/      generated outputs
data/                    ACIC 2022 challenge data, obtained separately (git ignored)
```

## Data

The ACIC 2022 Data Challenge datasets are distributed by the organisers at
<https://acic2022.mathematica.org/>. They are not redistributed here.

Unpack the Track 2 release so that the following paths exist, then either place it at
`./data` or point `ACIC22_DATA_ROOT` at it:

```
$ACIC22_DATA_ROOT/track2/practice/acic_practice_0001.csv       ... through 3400
$ACIC22_DATA_ROOT/track2/practice_year/acic_practice_year_0001.csv
$ACIC22_DATA_ROOT/ACIC_estimand_truths.csv
```

`ACIC_estimand_truths.csv` must contain the `dataset.num`, `variable`, `id.practice`
and `SATT` columns together with the four design factor columns, since the design
cells and the practice level truths are read from it.

## Environments

The two arms ran in two different environments, at different numpy major versions.
Install each in its own virtual environment and do not merge them.

```bash
python -m venv .venv-core
.venv-core/bin/pip install -r requirements-core.txt        # Scripts\ on Windows

python -m venv .venv-crosslib
.venv-crosslib/bin/pip install -r requirements-crosslib.txt
```

`requirements-core.txt` reproduces the environment recorded in `RUN_LOG.md` section 0,
which produced every confirmatory number in the manuscript: Python 3.11.0, numpy
1.26.4, pandas 2.3.3, scipy 1.15.3, scikit-learn 1.8.0, statsmodels 0.14.6, torch
2.11.0+cpu, dowhy 0.14, matplotlib 3.10.8.

`requirements-crosslib.txt` is a reconstruction rather than a bit exact recreation, and
the difference is worth stating plainly. The original cross library environment was a
plain virtual environment carrying DoubleML, EconML, catenets, jax and numpy 2.4.6,
with a path configuration file added by hand that appended the global site-packages
directory. numpy 2.4.6 therefore shadowed the global numpy 1.26.4 while scikit-learn
1.8.0, scipy 1.15.3, pandas 2.3.3, statsmodels 0.14.6 and torch 2.11.0+cpu were reused
from the global installation, where they had been built against numpy 1.26. Installing
the same version list cleanly from wheels, as this file does, instead yields binaries
built against the numpy 2 ABI. The cross library results should be treated as
reproducible to library version fidelity, not to binary fidelity.

That arrangement was not gratuitous. `RUN_LOG.md` section 124 records why: installing
into the global interpreter would have upgraded numpy or downgraded scikit-learn system
wide, a `--system-site-packages` virtual environment failed because a stale editable
install on `sys.path` broke pip's metadata scanner, and a fully isolated virtual
environment hit the Windows 260 character path limit while installing torch.

Reference hardware: Intel i7-1260P (16 threads), 16 GB RAM, no GPU, Windows 11. Wall
times below are from that machine and were affected by power throttling and by CPU
contention where the run log says so.

## Reproducing the analysis

Core environment, in order:

```bash
python code/validate_software.py                                          # checks V1 to V6
python code/run_fast_arm.py    --cohorts all --n-jobs 8 --chunk 200       # 3,400 cohorts, 48.1 min
python code/run_neural_arm.py  --per-cell 10 --seeds 5 --n-jobs 8 \
                               --chunk 10 --no-dowhy                      # 170 cohorts, 45.5 min
python code/run_e7_dowhy.py 10                                            # E7 cross check
python code/run_ovb_grid.py    --n-jobs 5 --adj reduced                   # 11x11 grid, 49.1 min
python code/analyze.py                                                    # tests and summaries
python code/assumption_checks_report.py
python code/make_report_tables.py
python code/figures.py                                                    # F1, F2, F4, F5, F6
python code/figures.py f3                                                 # F3, after the neural arm
python code/figure_fidelity_check.py                                      # figure QA
```

Cross library environment, after the core arm has written its tables:

```bash
python code/run_e3_e8_fast.py                                             # 3,400 cohorts x 2, 23.8 min
python code/run_catenets.py --seeds 5 --n-jobs 8 --chunk 16 --adj full    # 170 cohorts, 163.3 min
python code/analyze_addendum.py
```

Seeds are fixed throughout and no step depends on a global random state.
Cross fitting uses `StratifiedKFold(5, shuffle=True, random_state=42)` shared by E1, E2
and the placebo estimator; ridge `random_state` is 42; the cohort sample seed is
20260818; bootstrap seed 20260818 with B of 2,000, and 500 for concordance.
`run_neural_arm.py` and `run_catenets.py` checkpoint per chunk and resume.

## Software checks

V3 to V6 are self-contained and run from this repository alone. Recorded results, for
comparison against your own run:

- **V3** E-value at RR 3.9 equals 7.263, against the published worked example of 7.26.
- **V4** Rosenbaum `p_upper` at a bias parameter of 1 matches the scipy one sided
  normal approximation Wilcoxon p value to better than 1e-6, is monotone in the bias
  parameter, and is exact on a hand computed case.
- **V5** The bias bound reaches zero exactly at the robustness value. Linear robustness
  values match DoWhy 0.14 `linear-partial-R2` on cohorts 1 to 10, with maximum absolute
  differences of 2.45e-12 for RV, 2.44e-12 for RV at alpha, and 7.9e-10 for the estimate.
- **V6** On a synthetic linear design with n of 500 and a true effect of 2.0: TARNet
  plug-in 2.001, DragonNet plug-in 2.149, DragonNet doubly robust 2.039 ± 0.119.

**V1 and V2 need a second repository.** Both compare this code against a separate
reference implementation of the same two baseline estimators; V2 additionally reads that
implementation's `expected_results.json`. The reference implementation is the author's
benchmark repository, MIT licensed:

<https://github.com/harshsingh8005/rosenbound-benchmarks>

```bash
git clone https://github.com/harshsingh8005/rosenbound-benchmarks.git
export ACIC22_REPRODUCER_ROOT=/path/to/rosenbound-benchmarks/benchmarks/acic22_track2
export ACIC22_DATA_ROOT=/path/to/acic22            # both repositories read this variable
python code/validate_software.py
```

Left unset, V1 and V2 skip with a message and V3 to V6 still run. The relevant modules
are `src/aipw.py`, `src/dr_att.py` and `src/data_loader.py`; V2 reads
`expected_results.json` from the same directory.

Note that `src/rosenbaum.py` in that repository is a widening-bound surrogate,
`Γ* = |θ̂| / (z · se)`, and not a Rosenbaum randomisation inference bound. The frozen
plan records this and forbids describing it as one. The genuine matched pair Rosenbaum
bound used in this paper is in `code/sensitivity.py`, and the surrogate appears only as
a labelled comparator.

The manuscript reports both as passing: V1 with differences of exactly 0.0 in estimate
and standard error on cohorts 1, 500, 1000, 2000 and 3400, including identical covariate
matrices, column names, truths and trim counts; V2 reproducing the canonical 100 cohort
targets exactly. A reader without that implementation cannot independently re-verify
those two, and should treat V3 to V6 together with the independent DoubleML cross check
as the verifiable validation evidence. The DoubleML cross check is the stronger of the
two kinds in any case, since it compares against a third party library rather than
against another implementation by the same author.

## Deviations from the frozen plan

All four are recorded in `RUN_LOG.md` with timestamps. Summarised here because they
affect how the outputs should be read.

1. **Non-positive Riesz representer norm.** Plan section 5.2 defines the debiased form
   of the squared representer norm but did not cover the case where it is not positive,
   which occurs when the full covariate set induces near deterministic propensities.
   The adopted rule follows DoubleML's reference behaviour: fall back to the plug-in
   estimate and flag the cohort as `nu2_fallback`, with the fallback rate reported.
   Discovered during V5, before any confirmatory result, and so not outcome driven.
2. **Neural arm stall and logged re-run.** The pre-registered 90 minute cap expired with
   zero completed cohorts after a worker side stall at high parallelism. No neural result
   had been produced or inspected. The arm was re-run under the same pre-registered
   sample and seeds, and the DoWhy cross check was moved out of the parallel workers into
   its own single process script, which affects no pre-registered quantity.
3. **Figure F5 binning.** The pre-specified deciles of the product of the two measured
   strengths collapsed, because 6 of 10 deciles sat at exactly zero. F5 now uses one
   zero bin plus rank quintiles of the positive values, with an added paired difference
   panel. No hypothesis test depends on F5.
4. **Cross implementation replication outcome.** Under catenets reference defaults, with
   the same architectures, cohorts and seeds, the PyTorch neural conclusions do not
   replicate: the DragonNet coverage deficit is 5.3 points and not significant at
   p of 0.22, where the PyTorch run gave a confirmed 14.7 points, and the error ranking
   against the ridge baseline reverses. This is reported as the replication outcome, not
   corrected.

## Known limitation in the estimand

Every estimator in this pipeline weights each practice equally. There is no patient
weighting anywhere downstream of `acic_common.py`, while the challenge truths are the
scored sample average treatment effects on the treated. The manuscript names this as one
candidate mechanism for the coverage shortfall observed in the randomised arm, where no
confounding exists by construction. Anyone extending this code to a patient weighted
estimand should expect the randomised arm result to move.

## Citation

Please cite the manuscript rather than this repository. See `CITATION.cff` once the
preprint identifier is assigned.

## License

MIT. See `LICENSE`.
