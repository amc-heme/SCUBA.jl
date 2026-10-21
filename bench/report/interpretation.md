## Interpretation

These runs compare the same checksum-verified queries across SCUBA.jl, the SCUBA R package (legacy
R engine and Julia engine) and raw library baselines on synthetic datasets of 10k, 100k and 1M cells
(25k genes, ~1,500 nonzeros per cell). Numbers below are warm medians at 1M cells unless stated.

- **Storage layout decides feature-fetch cost.** Row-major (gene-compressed) BPCells makes a gene a
  contiguous read: SCUBA.jl fetches one gene from 1M cells in 2 ms and 1000 genes in 1.7 s. Col-major
  (cell-compressed) BPCells forces a full-matrix scan in every implementation, including BPCells' own
  C++ (3 s per query). csc h5ad behaves like row-major BPCells for single genes (5 ms) but scales worse
  with many genes (8.7 s for 1000); csr h5ad is unusable for feature queries at this scale (92 s).
  cellDIVER objects are row-major BPCells, so the row-major rows are the ones that matter for it.
- **On that layout the legacy R engine's cost is `SeuratObject::FetchData`, not storage.** Legacy
  Seurat + BPCells needs 2.3 s for one gene and 67 s for 1000; raw BPCells in R needs 85 ms and 4.7 s.
  The Julia engine brings the R package to 0.11 s and 4.7 s (20× and 14× faster than legacy) and
  sustains 30 queries/s across 8 Shiny-like sessions versus 5.4 for legacy.
- **Metadata-only calls go the other way through the Julia engine.** Tables round-trip R → Julia → R,
  so `fetch_metadata(full_table = TRUE)` costs 1.8 s versus 2 ms in legacy R, and `unique_values`,
  `fetch_cells` and `fetch_reduction` are 3–10× slower. cellDIVER's call mix is dominated by these,
  so the engine must route metadata and reduction calls to the R methods and send only expression
  fetches to Julia (planned next).
- **HDF5-backed SingleCellExperiment is a wash** between engines (0.2 s / 40 s): the slice is cut by
  DelayedArray in R either way, and 1000 genes over 1M cells is prohibitively slow on it.
- **Python.** Direct `h5py` slicing is the floor (0.1 ms per gene); backed `anndata` adds ~0.5 s of
  fixed overhead per query; the `scuba-sc` wrapper is bounded by building a pandas index of 1M cell
  IDs (158 ms) and by column conversion for wide results (12.9 s for 1000 genes), both planned fixes.
- **Concurrency.** HDF5's global lock caps in-process Julia threads at ~77 queries/s on h5ad while
  separate processes scale to 229 (csc) and 520 (row-major BPCells) queries/s at 16 sessions; BPCells
  has no such lock. For Shiny deployments this favours one engine process per session or per worker.
- **Cold start.** After evicting the page cache, a single-gene fetch costs 1.2–1.3 s on both h5ad and
  row-major BPCells (dominated by reading names and index structures), and opening a col-major
  BPCells source cold costs 8 s. Julia startup adds ~2 s and JuliaCall initialisation ~4 s before the
  first call; a sysimage and a precompile workload are the planned mitigations.
- **Memory.** A 1000-gene × 1M-cell dense result is 8 GB; peak RSS reached 25 GB in Julia and 60 GB
  in R. Streaming result columns, a `Float32` option and an enforceable size guard are planned.

Caveats: the legacy R anndata path was skipped on this machine (R `anndata` ≥ 0.8 returns objects
SCUBA 1.2.2 cannot dispatch on); csr h5ad was skipped beyond single genes at 1M; the 3M scale is
opt-in and was not run; timings are from one Apple M2 Ultra (24 cores, 206 GB RAM) with the datasets
on its internal SSD.
