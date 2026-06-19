# Simulation and Reproducibility Code

Run commands from the repository root.

The scripts in this folder support the current manuscript:

`Bayesian Model Averaging under Predictor Redundancy via Density-Ratio Posterior Compression`.

The commands below assume that `Rscript` is on the system path. If not, replace `Rscript` with the full path to the local `Rscript` executable.

## Quick Checks

```powershell
Rscript sim/run_support_kernel_unit_tests.R
Rscript sim/check_manuscript_artifacts.R
```

The unit test checks the density-ratio identities, safety constraint, pooled-pruned helper, active-set monotonicity, and related implementation diagnostics. The artifact check verifies that all table and figure inputs referenced by `tex/BMA.tex` exist.

## Main Benchmark Commands

Exact same-target and target-shift benchmark:

```powershell
Rscript sim/run_support_kernel_competitor_benchmark.R --mode=smoke
Rscript sim/run_support_kernel_competitor_benchmark.R --mode=full --checkpoint --workers=4 --chunk-size=8
Rscript sim/make_support_kernel_competitor_outputs.R --mode=full
```

Large non-enumerable support-space benchmark:

```powershell
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=smoke
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=full --checkpoint
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=summarize
Rscript sim/make_large_end_to_end_askpc_outputs.R
```

Reduced exact real-response spectroscopy benchmark:

```powershell
Rscript sim/run_real_response_reduced_benchmark.R --mode=smoke --dataset=all
Rscript sim/run_real_response_reduced_benchmark.R --mode=full --dataset=all
```

Split-stability and ablation checks:

```powershell
Rscript sim/run_split_stability.R --mode=screen
Rscript sim/run_adaptive_kernel_ablation.R --mode=medium
```

## Table and Figure Assembly

```powershell
Rscript sim/make_main_tables_figures.R
```

This assembles the current manuscript tables and figures from existing outputs. It does not rerun the expensive full benchmarks.

## Main Output Folders

- `sim/output/tables/` contains table inputs and their source CSVs.
- `sim/output/figures/` contains figure PDFs used by the manuscript.
- `sim/output/results/` contains compact source CSVs for large-\(p\) report diagnostics.
- `sim/output/large_end_to_end/` contains the retained reliable full-reference objects needed by the displayed large-\(p\) group-Hamming report.

Full benchmark runs are deterministic and checkpointed, but they can take many hours. Smoke runs are interface checks and should not be cited as completed manuscript evidence.
