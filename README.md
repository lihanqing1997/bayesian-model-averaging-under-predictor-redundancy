# Bayesian Model Averaging under Predictor Redundancy via Density-Ratio Posterior Compression

This companion repository contains the R code and output files used by the current paper.

The manuscript source is maintained separately. This repository is for reproducing the simulation tables, figures, diagnostics, and code-level checks.

The commands below assume that `Rscript` is on the system path. If not, replace `Rscript` with the full path to the local `Rscript` executable.

## Main Contents

- `sim/src/` contains support-kernel and density-ratio compression utilities.
- `sim/run_*.R` scripts run the exact, large-\(p\), reduced real-response, ablation, and stability benchmarks.
- `sim/make_*.R` and `sim/redraw_*.R` scripts assemble the current table and figure files.
- `sim/output/tables/` contains table inputs and source CSVs used by the current manuscript.
- `sim/output/figures/` contains figure PDFs used by the current manuscript.
- `sim/output/results/` contains compact source CSVs for large-\(p\) report diagnostics.
- `sim/output/large_end_to_end/` contains retained reliable full-reference objects needed by the displayed large-\(p\) group-Hamming report.

## Quick Checks

```bash
Rscript sim/run_support_kernel_unit_tests.R
```

The artifact checker `sim/check_manuscript_artifacts.R` verifies table and figure inputs against a manuscript TeX file. Use it from the manuscript-source checkout, or pass the TeX file path explicitly as the first argument.

## Main Benchmarks

Exact same-target and target-shift benchmark:

```bash
Rscript sim/run_support_kernel_competitor_benchmark.R --mode=smoke
Rscript sim/run_support_kernel_competitor_benchmark.R --mode=full --checkpoint --workers=4 --chunk-size=8
Rscript sim/make_support_kernel_competitor_outputs.R --mode=full
```

Large non-enumerable support-space benchmark:

```bash
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=smoke
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=full --checkpoint
Rscript sim/run_large_end_to_end_askpc_benchmark.R --mode=summarize
Rscript sim/make_large_end_to_end_askpc_outputs.R
```

Reduced exact real-response spectroscopy benchmark:

```bash
Rscript sim/run_real_response_reduced_benchmark.R --mode=smoke --dataset=all
Rscript sim/run_real_response_reduced_benchmark.R --mode=full --dataset=all
```

Table and figure assembly:

```bash
Rscript sim/make_main_tables_figures.R
```

## Reading the Diagnostics

- `TV`, `FKL`, and `RKL` measure posterior distortion relative to the chosen reference support posterior.
- `q0` is the fallback weight on the unrestricted posterior. Small distortion with large `q0` is only partial compression.
- `expected_code` is the internal column name for expected reporting cost.
- `active_kernels_001`, `q_effective_kernels`, `n_kernels`, and `stored_atoms` are list-level reporting diagnostics.

Dilution-prior and DPP-prior rows change the Bayesian target. Their TV/FKL values should be read as target-shift diagnostics rather than same-posterior compression failures.
