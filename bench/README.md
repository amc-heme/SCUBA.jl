# SCUBA benchmark suite

Reproducible performance comparison of SCUBA's frontends on disk-backed single-cell storage:

- **Julia** `SCUBA.jl` native readers: backed `.h5ad` (csc and csr), BPCells matrix directories, in-memory snapshot
- **R** `SCUBA` package: legacy R engine and the Julia engine (`options(SCUBA.engine = "julia")`) on Seurat v5 objects with BPCells layers and on HDF5-backed SingleCellExperiment objects
- **Baselines** without SCUBA: BPCells in R, `anndata` backed mode and direct `h5py` slicing in Python, the `scuba-sc` Python package

on synthetic datasets of 10k, 100k, 1M and (opt-in) 3M cells × 25k genes with ~1,500 nonzeros per cell.
Every timed call is checksummed against values computed by the generator, so timings are only ever
compared across implementations that returned identical data.

## Layout

| path | purpose |
|---|---|
| `config/bench.toml` | scales, workloads, repetitions, skip rules, concurrency mix |
| `generate/make_dataset.py` | synthetic data → `csc.h5ad`, `csr.h5ad`, parquet tables, `selection.json`, `expected.json` |
| `generate/make_r_objects.R` | BPCells directories, Seurat RDS, HDF5-backed SCE RDS |
| `julia/run.jl`, `r/run.R`, `python/run.py` | harnesses (one process per language/implementation/scale) |
| `julia/worker.jl`, `r/worker.R`, `python/worker.py`, `tools/concurrency.py` | concurrent-session benchmark |
| `tools/cold.sh`, `tools/run_cold.py` | page-cache eviction (`vmtouch -e`) and cold runs in fresh processes |
| `tools/preflight.py`, `tools/run_with_timeout.py`, `tools/merge_results.py` | checks, timeouts, merging raw CSVs to Parquet |
| `report/make_report.R`, `report/report.Rmd` | figures and `summary.md` from `results/results.parquet` |
| `common/schema.md` | the result-row contract |

Datasets live in `~/scuba-bench-data/<scale>/` (`[suite].data_dir`), results in `results/raw/<run_id>/`.

## Requirements

- Julia ≥ 1.10 (this project: `julia --project=bench -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`)
- conda env `scuba-fixtures` (`conda env update -n scuba-fixtures -f bench/environment.yml`)
- R with SeuratObject, Seurat, BPCells, SingleCellExperiment, HDF5Array, JuliaCall, arrow, ps, RcppTOML, digest, ggplot2, dplyr, tidyr, rmarkdown; the SCUBA R package checkout at `[suite].r_scuba_path` (branch `feature/julia-engine`)
- `vmtouch` (`brew install vmtouch`) for cold-cache runs

## Running

```sh
cd bench
make preflight SCALE=1M          # disk / RAM / tools; writes results/env/<run_id>.json
make smoke                       # 10k: generate, Julia + Python harness, merge, report (R=1 adds the R harness)
make data-100k r-objects-100k    # build a dataset and its R objects
make julia-100k r-100k python-100k RUN_ID=my-run
make cold-100k conc-100k RUN_ID=my-run
make report                      # results/results.parquet -> report/figures, report/summary.md, report/report.html
make all-1M RUN_ID=my-run        # everything at 1M (hours)
make data-3M                     # opt-in; ~180 GB on disk, ~120 GB RAM to generate
```

Harness flags: `--only impl1,impl2`, `--resume` (skip measured cells), `--include-optional` (in-memory
baselines), `--single "feature|data|10:0|all" --rep N` (one workload once, used by the cold driver),
`--sweep max_read_elements=1048576,8388608` (Julia read-size sweep).

## Method

- Warm: one discarded first call (recorded as `phase = first_call`, it includes JIT) then 5 repetitions
  per cell; cells whose first call exceeds `reduce_reps_above_s` get a single repetition; cells predicted
  (from the previous scale) to exceed `timeout_s` are skipped as `skipped_extrapolated`.
- Cold: `vmtouch -e` over the storage's files, verification with `vmtouch -v`, one fresh process per
  repetition; unverified evictions are labelled `cold_unverified`.
- Concurrency: N sessions (processes, or Julia threads) run a weighted query mix for `duration_s`
  after a barrier; throughput and p50/p95/p99 per workload go to `concurrency_summary.csv`.
- Memory: Julia `@timed` bytes and `Sys.maxrss`; R `gc()` max-used and `ps` RSS; Python `tracemalloc`
  (≤ 100k) and `getrusage`.
- All sources are read with the same gene sets (`selection.json`, seeded, identical across scales) and
  cell subsets (`all`, 1% / 10% random, 1% / 10% contiguous).
