## Interpretation

These runs compare the same checksum-verified queries across SCUBA.jl, the SCUBA R package (legacy
R engine and Julia engine) and raw library baselines on synthetic datasets of 10k, 100k and 1M cells
(25k genes, ~1,500 nonzeros per cell). Numbers below are warm medians at 1M cells unless stated.

- **Storage layout decides feature-fetch cost.** Row-major (gene-compressed) BPCells makes a gene a
  contiguous read: SCUBA.jl fetches one gene from 1M cells in 2 ms and 1000 genes in 1.7 s. Col-major
  (cell-compressed) BPCells forces a full-matrix scan in every implementation, including BPCells' own
  C++ (3 s per query). csc h5ad is the fastest layout for feature fetches once the reader is
  type-stable (1.2 ms for one gene, 1.0 s for 1000; an earlier version of the reader boxed every
  nonzero and needed 5 ms / 8.7 s); csr h5ad still scans the whole matrix (3.7 s per query, skipped
  beyond single genes). cellDIVER objects are row-major BPCells, so the row-major rows are the ones
  that matter for it.
- **On that layout the legacy R engine's cost is `SeuratObject::FetchData`, not storage.** Legacy
  Seurat + BPCells needs 2.3 s for one gene and 67 s for 1000; raw BPCells in R needs 85 ms and 4.7 s.
  The Julia engine brings the R package to 0.11 s and 4.7 s (20× and 14× faster than legacy) and
  sustains 30 queries/s across 8 Shiny-like sessions versus 5.4 for legacy.
- **Metadata-only calls went the other way through the Julia engine.** Tables round-trip R → Julia → R,
  so `fetch_metadata(full_table = TRUE)` cost 1.8 s versus 2 ms in legacy R, and `unique_values`,
  `fetch_cells` and `fetch_reduction` were 3–10× slower. cellDIVER's call mix is dominated by these,
  so the R package now routes metadata and reduction calls to the R methods and sends only
  expression fetches to Julia; the `[julia]` rows for those workloads therefore equal the legacy
  rows in current builds (the tables above still show the pre-routing measurements).
- **HDF5-backed SingleCellExperiment is a wash** between engines (0.2 s / 40 s): the slice is cut by
  DelayedArray in R either way, and 1000 genes over 1M cells is prohibitively slow on it.
- **Python.** Direct `h5py` slicing is the floor (0.1 ms per gene, 128 ms for 1000); backed `anndata`
  adds ~0.5 s of fixed overhead per query; the `scuba-sc` wrapper costs 7.7 ms for one gene and 3.3 s
  for 1000 (it was 158 ms and 12.9 s before the cell index was cached as an Arrow-backed index and
  wide results were handed to pandas as one zero-copy block). The remaining gap to `h5py` is the
  pandas frame itself: an 8 GB block has to be gathered once.
- **Concurrency.** HDF5's global lock caps in-process Julia threads at ~77 queries/s on h5ad while
  separate processes scale to 229 (csc) and 520 (row-major BPCells) queries/s at 16 sessions; BPCells
  has no such lock. For Shiny deployments this favours one engine process per session or per worker.
- **Cold start.** After evicting the page cache, a single-gene fetch costs 1.2–1.3 s on both h5ad and
  row-major BPCells (dominated by reading names and index structures), and opening a col-major
  BPCells source cold costs 8 s. Process start-up (measured separately, 100k objects): Julia 1.4 s →
  0.26 s with a PackageCompiler sysimage; Python `scuba-sc` 5.9 s → 2.1 s; the R engine's Julia
  start-up 5.0 s → 3.3 s (JuliaCall's own set-up is ~2 s of that and RCall cannot be compiled into
  the image on Julia 1.13), after which the first call on a real object costs ~0.15 s thanks to a
  warm-up on tiny objects at start-up. For a Shiny app this is a one-time cost in `global.R`.
- **File-backed h5ad from R (`r_h5ad_julia`).** `SCUBA::open_h5ad()` serves a csc h5ad file to R
  through the Julia engine without reticulate: one gene over 1M cells costs 245 ms (13 ms at 100k)
  and 1000 genes 4.1 s; `fetch_matrix(sparse = TRUE)` returns the 1000-gene block as a `dgCMatrix`
  in 3.1 s. The fixed cost per call is the 1M-element cell-name column crossing JuliaCall into R
  (the Julia reader itself needs 1.2 ms); metadata, cells and reductions stay in R (81 ms for the
  full table, 87 ms for `fetch_cells`). This replaces the legacy reticulate path, which could not
  be measured at all with current R `anndata`.
- **Sparse and chunked results in Python.** `fetch_matrix(sparse=True)` hands 1000 genes × 1M cells
  to scipy in 1.0 s versus 3.1 s for the dense pandas frame; `iter_chunks` over the same request
  (sparse, 256 MiB chunks) takes 10 s because every chunk re-reads the requested csc columns, so
  chunking is for bounding memory, not for speed, on this layout.
- **R bridge overheads removed (2026-10-22, evening).** Cell names no longer cross the JuliaCall
  bridge in either direction (R already holds them), and expression columns are written by Julia
  straight into R double vectors instead of being built in Julia and copied. One gene over 1M
  cells from R: file-backed h5ad 245 ms → 26 ms, row-major Seurat 110 ms → 16 ms (bench medians).
  1000 genes: h5ad 4.1 s → 2.9 s and Seurat 4.7 s → 2.7 s in the harness, where every repetition
  allocates and touches 8 GB of fresh R memory; a warm R session that reuses that memory does the
  h5ad fetch in 0.9 s. Process peak for the 1000-gene fetch: Seurat 18.6 GB → 12.0 GB, h5ad
  10.6 GB, against the 8 GB result itself. `options(SCUBA.julia_eltype)` is therefore a no-op.
- **Memory.** A 1000-gene × 1M-cell dense result is 8 GB; peak RSS reached 25 GB in Julia and 60 GB
  in R in these runs. SCUBA.jl now fills result columns directly for results above 256 MiB (peak
  ~9 GB for the same query), offers `eltype = Float32` and a sparse `fetch_matrix`, and the size
  guard can be switched from a warning to an error for deployments.

Caveats: the legacy R anndata path was skipped on this machine (R `anndata` ≥ 0.8 returns objects
SCUBA 1.2.2 cannot dispatch on); csr h5ad was skipped beyond single genes at 1M; the 3M scale is
opt-in and was not run; timings are from one Apple M2 Ultra (24 cores, 206 GB RAM) with the datasets
on its internal SSD.
