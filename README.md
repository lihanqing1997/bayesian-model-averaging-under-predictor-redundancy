# Bayesian Model Averaging under Predictor Redundancy

This repository contains the simulation and reproducibility code for
`Bayesian Model Averaging under Predictor Redundancy`.

The manuscript LaTeX source, bibliography, and journal style files are not
included in this public companion repository. During review, the manuscript is
maintained separately. This repository is intended to reproduce the empirical
tables, figures, diagnostics, and code-level checks used by the paper.

The submitted paper is centered on density-ratio support-kernel posterior
compression for a fixed unrestricted support posterior. The theory is written
for support-indexed Bayesian regression models in general. The reproducible
experiments use a regularized Gaussian linear BMA reference posterior because it
supports exact enumeration, reference-sampler diagnostics, and direct
competitor comparisons. Older representative-family, direct sampler,
retained-mass, and discarded empirical material is kept only as archived backup
or historical generated output. It is not part of the active paper.

## Repository Contents

- `sim/src/support_kernels.R` contains support-kernel utilities.
- `sim/src/adaptive_support_kernel_compression.R` contains active-set compression routines.
- `sim/output/tables/` contains table source files.
- `sim/output/figures/` contains figure source files.
- `sim/run_*.R` and `sim/make_*.R` scripts regenerate benchmark outputs and
  assembled table or figure files.

## How To Read Compression Results

The paper treats density-ratio support-kernel posterior compression as a
reporting layer for a fixed unrestricted support posterior. It is not a new
prior, not a universal variable-selection rule, and not a replacement for the
reference posterior sampler.

The main compression diagnostics must be read together.

- `TV`, `FKL`, and `RKL` measure posterior distortion relative to the fixed
  reference posterior.
- `q0` is the fallback weight on the unrestricted posterior. Small distortion
  with large `q0` is only partial compression.
- `expected_code` is the expected descriptor cost
  `sum_m q_m c_m` used by the optimization objective.
- `active_kernels_001`, `q_effective_kernels`, `n_kernels`, and `stored_atoms`
  are list-level storage checks. The storage-sensitivity table also reports a
  `fallback_reading` column derived from `q0`.
- In the large `p=100` table, `Storage` is expected descriptor cost and `List`
  is the active-kernel or stored-atom count used as a compact list-level check.

Dilution-prior BMA and DPP-prior BMA change the posterior target. In exact
enumerable experiments their TV and FKL columns are target-shift diagnostics.
The dilution row uses a powered determinant correction and the DPP row uses a
stronger determinant prior weight. In the large `p=100` runs, those rows are
reweighted over the support set visited by the unrestricted reference sampler.
A production changed-prior benchmark would need a separate changed-prior
sampler or importance-weight effective-sample-size diagnostics.

## R Requirements

The scripts were developed with R 4.4.2. Install the main packages with

```r
install.packages(
  c("yaml", "glmnet", "grpreg", "pls", "coda", "boot", "SGL", "ggplot2"),
  repos = "https://cloud.r-project.org"
)
```

The commands below assume that `Rscript` is available on the system path. If it
is not, replace `Rscript` with the full path to the local Rscript executable.

## Main Support-Kernel Runs

Run the fast code-level checks first. They verify exact density-ratio
identities, the safety constraint, active-set monotonicity, oracle-cover
diagnostics, blocked Monte Carlo summaries, and the pooled-pruned support-kernel
helper.

```bash
Rscript sim/run_support_kernel_unit_tests.R
```

Smoke mode checks interfaces. Medium mode is a development benchmark. Full mode
is deterministic, resumable, and regenerates the current exact competitor table
and rate-distortion figure used in Section 6.

```bash
Rscript sim/run_main_support_kernel_smoke.R
Rscript sim/run_main_support_kernel_medium.R
Rscript sim/run_main_support_kernel_full.R
```

The next benchmark layer compares support-kernel dictionary constructions with
the closest empirical comparison families from the revised Section 6. It uses
unrestricted exact BMA as the reference target and evaluates Top-M support
atoms, credible support sets, posterior-cluster kernels, dilution-prior BMA,
DPP-prior BMA, fixed hard dictionaries, and pooled-pruned support-kernel
dictionaries under the same diagnostic schema. Lasso, elastic net, and group
lasso are included as prediction-only baselines. They are not assigned TV or FKL
because they are not posterior distributions over supports.

```bash
Rscript sim/run_support_kernel_competitor_benchmark.R --mode=smoke
Rscript sim/run_support_kernel_competitor_benchmark.R --mode=pilot
Rscript sim/run_support_kernel_competitor_benchmark.R --mode=medium
Rscript sim/run_support_kernel_competitor_benchmark.R --mode=full --checkpoint --workers=4 --chunk-size=8
Rscript sim/make_support_kernel_competitor_outputs.R --mode=full
Rscript sim/run_logistic_support_kernel_check.R --mode=full
Rscript sim/make_q0_screened_frontier_outputs.R
```

The exact benchmark and the large end-to-end benchmark answer different
questions. The exact benchmark is synthetic and small enough to enumerate all
supports, so TV/FKL/RKL are free of MCMC error. The large end-to-end benchmark is
also synthetic, but non-enumerable. It first generates data, runs multi-chain
support MCMC to obtain an empirical unrestricted BMA reference posterior,
records reference diagnostics, and then compares pooled-pruned support kernels
with Top-M atoms, credible support sets, posterior-cluster kernels,
dilution-prior BMA and DPP-prior BMA over the visited support set, and fixed
hard dictionaries. Lasso, elastic net, and group lasso are again reported only
as predictive baselines.
Smoke and pilot modes are interface checks. Full mode now uses a longer
six-chain reference run with 20,000 iterations, 5,000 burn-in iterations, and
5-step thinning by default. It is checkpointed and can be run in shards.

```bash
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=smoke
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=pilot
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=full --checkpoint
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=summarize --publish
Rscript sim/make_large_end_to_end_askpc_outputs.R
```

The reported large-benchmark tables use the `reliable_full_v2` full MCMC shard outputs
augmented with credible-set, DPP-prior, lasso, elastic-net, and group-lasso
competitor rows. To rebuild those augmented rows from the saved reference
objects, run

```bash
Rscript sim/augment_large_end_to_end_competitors.R --run-tag=reliable_full_v2 --out-tag=reliable_full_v2_augmented --beta-grid=0.02 --coverage-grid=0.9 --dpp-power-grid=0.5,1
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=summarize --checkpoint-path=sim/output/tables/table_large_end_to_end_askpc_detail_reliable_full_v2_augmented.csv --publish
Rscript sim/make_large_end_to_end_askpc_outputs.R
```

Useful full-mode overrides are `--replications=5`, `--scenarios=one_representative,weak_signal`,
`--rho=0.7,0.9`, `--n-iter=20000`, `--chains=6`, `--shard-id=1`, and
`--shard-total=4`. The shard options allow several independent R processes to
work on different benchmark cells. Add `--run-tag=name` to keep a new rerun
separate from the released result tables while it is being checked. Summarize the
tagged checkpoints with `--mode=summarize --run-tag=name`, and add `--publish`
only after the diagnostic table has been inspected.

The adaptive-kernel runners are also available as convenience wrappers.

```bash
Rscript sim/run_all_adaptive_kernel_smoke.R
Rscript sim/run_all_adaptive_kernel_medium.R
Rscript sim/run_all_adaptive_kernel_full.R
```

The reduced exact real-response spectroscopy benchmark uses real Tecator and
Gasoline responses and real spectral designs. It averages contiguous spectral
channels into 12 ordered regions in full mode using an X-only rule and then
enumerates the exact BMA posterior. This is real-response evidence under an
exact reference posterior, not a full-resolution MCMC stress test. In this
benchmark the dilution-prior row uses the powered determinant correction and
the DPP-prior row uses unit determinant power, so the two changed-prior
comparisons are intentionally distinct.

```bash
Rscript sim/run_real_response_reduced_benchmark.R --mode=smoke --dataset=all
Rscript sim/run_real_response_reduced_benchmark.R --mode=full --dataset=all
```

## Table And Figure Assembly

Regenerate the current main tables and figures with

```bash
Rscript sim/run_adaptive_kernel_ablation.R --mode=medium
Rscript sim/make_storage_code_outputs.R
Rscript sim/make_minor_revision_outputs.R
Rscript sim/make_semisynthetic_tecator_kernel_figure.R
Rscript sim/redraw_main_figures.R
Rscript sim/make_main_tables_figures.R
```

## Generated Table And Figure Outputs

The reported analysis uses the following table outputs.

- `sim/output/tables/table_soft_kernel_identity.csv`
- `sim/output/tables/table_topm_atom_separation_summary.csv`
- `sim/output/tables/table_adaptive_kernel_exact.csv`
- `sim/output/tables/table_adaptive_kernel_ablation.csv`
- `sim/output/tables/table_adaptive_kernel_ablation_pooled_pruned.csv`
- `sim/output/tables/table_adaptive_kernel_ablation_convergence.csv`
- `sim/output/tables/table_adaptive_kernel_full_draws.csv`
- `sim/output/tables/table_adaptive_kernel_semisynthetic_realx.csv`
- `sim/output/tables/table_semisynthetic_tecator_band_overlap.csv`
- `sim/output/tables/table_support_kernel_competitor_benchmark_detail.csv`
- `sim/output/tables/table_support_kernel_competitor_benchmark_summary.csv`
- `sim/output/tables/table_support_kernel_competitor_best_fkl.csv`
- `sim/output/tables/table_support_kernel_algorithm_diagnostics.csv`
- `sim/output/tables/table_support_kernel_predictive_baselines.csv`
- `sim/output/tables/table_logistic_support_kernel_check.csv`
- `sim/output/tables/table_logistic_support_kernel_check_detail.csv`
- `sim/output/tables/table_q0_screened_frontier.csv`
- `sim/output/tables/table_storage_code_definitions.csv`
- `sim/output/tables/table_storage_code_sensitivity.csv`
- `sim/output/tables/table_large_end_to_end_askpc_detail.csv`
- `sim/output/tables/table_large_end_to_end_askpc_summary.csv`
- `sim/output/tables/table_large_end_to_end_askpc_best_fkl.csv`
- `sim/output/tables/table_large_end_to_end_reference_diagnostics.csv`
- `sim/output/tables/table_large_end_to_end_askpc_family_summary.csv`
- `sim/output/tables/table_large_end_to_end_reference_diagnostics_summary.csv`
- `sim/output/tables/table_large_end_to_end_predictive_baselines.csv`
- `sim/output/tables/table_key_paired_comparisons.csv`
- `sim/output/tables/table_reference_computation_cost.csv`
- `sim/output/tables/table_real_response_reduced_detail.csv`
- `sim/output/tables/table_real_response_reduced_summary.csv`
- `sim/output/tables/table_real_response_reduced_best_fkl.csv`
- `sim/output/tables/table_real_response_reduced_reference_diagnostics.csv`
- `sim/output/figures/fig_support_kernel_competitor_frontier.pdf`

The reported analysis uses the following figure outputs.

- `sim/output/figures/fig_topm_atom_separation.pdf`
- `sim/output/figures/fig_adaptive_kernel_rate_distortion.pdf`
- `sim/output/figures/fig_adaptive_kernel_semisynthetic_realx.pdf`
- `sim/output/figures/fig_semisynthetic_tecator_kernels.pdf`

## Theory And Diagnostic Utilities

The code package includes utilities that correspond to the main
validation results.

- `kernel_pool_cover_diagnostics()` in `sim/src/support_kernels.R` computes the
  normalized-column cover distance used by the oracle-cover theorem.
- `kernel_code_length()`, `near_duplicate_kernel_keep()`, and
  `kernel_summary_table()` in `sim/src/support_kernels.R` define the common
  code-length interface, near-duplicate normalized-column pruning, and the top
  selected-kernel summaries reported by the pooled-pruned support-kernel fit.
- `sim/make_storage_code_outputs.R` writes the reporting-code convention table
  and a storage-sensitivity summary that contrasts expected descriptor code
  with active-list, total-list, and stored-atom diagnostics where available.
- `posterior_functional_error()` in
  `sim/src/family_mixture_compression.R` checks the TV bound for bounded
  posterior functionals.
- `blocked_validation_summary()` in
  `sim/src/family_mixture_compression.R` reports block-based Monte Carlo
  standard errors for Markov-chain validation output.
- `fit_askpc_pooled_pruned()` in
  `sim/src/adaptive_support_kernel_compression.R` implements the default
  pooled-pruned fitting path. It fits the pooled dictionary, prunes by
  non-safety q-mass, refits, warm-starts the pooled problem from the pruned
  solution, and checks the pooled-versus-pruned objective relation to numerical
  tolerance.
- `simplex_directional_kkt_residual()` in
  `sim/src/family_mixture_compression.R` reports the feasible-direction KKT
  residual used by the optimizer diagnostics.

## Evidence Labels

The paper distinguishes the following evidence types.

- Exact enumerable checks verify the density-ratio identities.
- Exact simulations use synthetic data with small enough \(p\) to enumerate all
  \(2^p\) supports, so posterior-compression metrics have no MCMC error.
- Medium-mode validations test the adaptive dictionary and ablation workflow.
- Saved-posterior-draw compression evaluates the method on stored unrestricted BMA draws.
- Large end-to-end simulations use non-enumerable synthetic data, run the
  reference-posterior sampler, report diagnostics, and then run compression.
  They must be read together with the reported reference diagnostics.
  The active full \(p=100\) rerun uses the `reliable_full_v2` shard outputs and passes the core reference-posterior gate in all 36 cells, with one low-acceptance warning.
- Semi-synthetic real-\(X\) validation preserves predictor geometry while using a known simulated signal.
- Reduced exact real-response spectroscopy uses real predictors and real responses,
  with exact enumerated reference posteriors after deterministic spectral
  region averaging.

Do not use smoke or debug outputs as completed manuscript evidence.

## Random Seeds

The global seed base is in `sim/config_jmlr.yaml`. Seeds are deterministic
functions of scenario, correlation, replication, method, and run mode.

## Archived Material

Older material from previous versions is preserved under `archive/` and in
historical output folders. Those files are useful for provenance, but they are
not included in the reported results and should not be cited as current
evidence.
