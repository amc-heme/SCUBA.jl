# SCUBA benchmark summary

Generated 2026-10-20 23:34 MDT from `results.parquet` (844719 rows, 843976 ok).

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

## 10k cells: feature fetch over all cells, warm (median, IQR)

| implementation | genes=1 | genes=10 | genes=100 | genes=1000 |
|---|---|---|---|---|
| jl_h5ad_csc | 0.2 ms (IQR 0.1 ms–0.2 ms) | 1.1 ms (IQR 1.0 ms–1.7 ms) | 11.7 ms (IQR 11.1 ms–12.6 ms) | 125.5 ms (IQR 124.0 ms–126.0 ms) |
| jl_h5ad_csr | 890.7 ms (IQR 880.2 ms–904.4 ms) | 896.2 ms (IQR 893.2 ms–908.5 ms) | 980.4 ms (IQR 966.6 ms–997.4 ms) | 1.0 s (IQR 987.6 ms–1.0 s) |
| jl_bpcells | 30.5 ms (IQR 30.4 ms–30.6 ms) | 30.8 ms (IQR 30.7 ms–30.9 ms) | 32.7 ms (IQR 32.6 ms–33.0 ms) | 42.3 ms (IQR 42.2 ms–42.6 ms) |
| jl_bpcells_rowmajor | 0.0 ms (IQR 0.0 ms–0.0 ms) | 0.1 ms (IQR 0.1 ms–0.1 ms) | 1.0 ms (IQR 0.9 ms–1.0 ms) | 10.3 ms (IQR 9.9 ms–10.3 ms) |
| raw_anndata_py | 8.6 ms (IQR 8.6 ms–8.8 ms) | 8.8 ms (IQR 8.5 ms–8.9 ms) | 13.3 ms (IQR 12.8 ms–13.7 ms) | 44.8 ms (IQR 44.3 ms–45.6 ms) |
| raw_h5py_py | 0.1 ms (IQR 0.1 ms–0.1 ms) | 1.2 ms (IQR 1.1 ms–1.2 ms) | 9.9 ms (IQR 9.9 ms–10.5 ms) | 100.4 ms (IQR 99.8 ms–105.7 ms) |
| py_scuba_sc | 2.7 ms (IQR 2.6 ms–2.7 ms) | 4.1 ms (IQR 3.8 ms–4.7 ms) | 24.9 ms (IQR 22.9 ms–33.5 ms) | 254.5 ms (IQR 244.5 ms–267.2 ms) |
| r_seurat_bpcells [r] | 54.0 ms (IQR 52.5 ms–55.0 ms) | 61.0 ms (IQR 60.0 ms–62.0 ms) | 125.0 ms (IQR 123.5 ms–130.0 ms) | 1.3 s (IQR 1.3 s–1.3 s) |
| r_seurat_bpcells_rowmajor [r] | 26.0 ms (IQR 25.5 ms–26.0 ms) | 33.0 ms (IQR 32.0 ms–33.0 ms) | 102.0 ms (IQR 99.5 ms–103.5 ms) | 1.4 s (IQR 1.3 s–1.4 s) |
| r_sce_backed [r] | 35.0 ms (IQR 35.0 ms–37.5 ms) | 38.0 ms (IQR 36.0 ms–38.0 ms) | 57.0 ms (IQR 56.0 ms–57.0 ms) | 542.0 ms (IQR 248.5 ms–571.0 ms) |
| raw_bpcells_r [r] | 35.0 ms (IQR 34.0 ms–36.5 ms) | 35.0 ms (IQR 34.0 ms–36.0 ms) | 44.0 ms (IQR 42.0 ms–46.0 ms) | 80.0 ms (IQR 79.0 ms–82.0 ms) |
| raw_bpcells_r_rowmajor [r] | 3.0 ms (IQR 3.0 ms–4.0 ms) | 5.0 ms (IQR 4.0 ms–5.0 ms) | 13.0 ms (IQR 12.0 ms–15.5 ms) | 91.0 ms (IQR 89.0 ms–99.0 ms) |
| r_seurat_bpcells [julia] | 32.0 ms (IQR 30.5 ms–32.0 ms) | 32.0 ms (IQR 32.0 ms–33.0 ms) | 39.0 ms (IQR 38.0 ms–41.0 ms) | 127.0 ms (IQR 110.5 ms–135.0 ms) |
| r_seurat_bpcells_rowmajor [julia] | 1.0 ms (IQR 1.0 ms–1.0 ms) | 2.0 ms (IQR 1.0 ms–2.0 ms) | 7.0 ms (IQR 6.0 ms–7.0 ms) | 102.0 ms (IQR 79.5 ms–107.5 ms) |
| r_sce_backed [julia] | 9.0 ms (IQR 9.0 ms–9.0 ms) | 12.0 ms (IQR 12.0 ms–13.0 ms) | 48.0 ms (IQR 45.0 ms–49.0 ms) | 500.0 ms (IQR 481.0 ms–764.5 ms) |

## 10k cells: metadata, reductions, cells (median)

| implementation | fetch_cells:batch:b3 | fetch_cells:cell_type:T_cell | metadata:cell_type | metadata:full | metadata:nCount_RNA | reduction:X_pca|1-2 | reduction:X_pca|1-50 | reduction:X_umap|1-2 | unique_values:cell_type | unique_values:sample |
|---|---|---|---|---|---|---|---|---|---|---|
| jl_h5ad_csc | 1.5 ms | 1.6 ms | 0.7 ms | 4.4 ms | 0.1 ms | 2.0 ms | 44.3 ms | 1.8 ms | 0.3 ms | 0.3 ms |
| jl_h5ad_csr | 1.4 ms | 1.5 ms | 0.7 ms | 4.4 ms | 0.1 ms | 2.0 ms | 44.6 ms | 1.8 ms | 0.3 ms | 0.3 ms |
| jl_bpcells | 1.1 ms | 1.2 ms | 0.4 ms | 0.2 ms | 0.0 ms | 0.0 ms | 0.4 ms | 0.0 ms | 0.0 ms | 0.0 ms |
| jl_bpcells_rowmajor | 1.4 ms | 1.5 ms | 0.4 ms | 0.2 ms | 0.0 ms | 0.0 ms | 0.4 ms | 0.0 ms | 0.0 ms | 0.0 ms |
| raw_anndata_py | 9.4 ms | 9.8 ms | 19.4 ms | 0.1 ms | 0.1 ms | 0.1 ms | 0.4 ms | 0.0 ms | 19.9 ms | 19.1 ms |
| py_scuba_sc | 2.5 ms | 3.6 ms | 23.7 ms | 12.3 ms | 3.3 ms | 5.2 ms | 52.7 ms | 4.9 ms | 0.4 ms | 0.4 ms |
| r_seurat_bpcells [r] | 1.0 ms | 2.0 ms | 2.0 ms | 0.0 ms | 1.0 ms | 1.0 ms |  | 0.0 ms | 1.0 ms | 1.0 ms |
| r_seurat_bpcells_rowmajor [r] | 1.0 ms | 1.0 ms | 2.0 ms | 0.0 ms | 1.0 ms | 1.0 ms |  | 1.0 ms | 1.0 ms | 1.0 ms |
| r_sce_backed [r] | 2.0 ms | 2.0 ms | 4.0 ms | 6.0 ms | 3.0 ms | 19.0 ms |  | 19.0 ms | 2.0 ms | 3.0 ms |
| r_seurat_bpcells [julia] | 2.0 ms | 3.0 ms | 4.0 ms | 14.0 ms | 2.0 ms | 2.0 ms | 14.0 ms | 2.0 ms | 4.0 ms | 4.0 ms |
| r_seurat_bpcells_rowmajor [julia] | 3.0 ms | 3.0 ms | 4.0 ms | 16.0 ms | 1.0 ms | 2.0 ms | 15.0 ms | 2.0 ms | 4.0 ms | 4.0 ms |
| r_sce_backed [julia] | 2.0 ms | 3.0 ms | 3.0 ms | 19.0 ms | 2.0 ms | 5.0 ms | 18.0 ms | 5.0 ms | 4.0 ms | 4.0 ms |

## 100k cells: feature fetch over all cells, warm (median, IQR)

| implementation | genes=1 | genes=10 | genes=100 | genes=1000 |
|---|---|---|---|---|
| jl_h5ad_csc | 0.6 ms (IQR 0.5 ms–4.8 ms) | 7.3 ms (IQR 4.7 ms–7.9 ms) | 83.9 ms (IQR 77.2 ms–98.9 ms) | 909.5 ms (IQR 887.2 ms–939.8 ms) |
| jl_h5ad_csr | 9.9 s (IQR 9.8 s–9.9 s) | 9.9 s (IQR 9.9 s–10.0 s) | 9.5 s (IQR 9.4 s–9.5 s) | 9.9 s (IQR 9.8 s–10.0 s) |
| jl_bpcells | 317.3 ms (IQR 308.2 ms–322.9 ms) | 315.9 ms (IQR 305.4 ms–324.0 ms) | 367.6 ms (IQR 365.2 ms–371.8 ms) | 407.6 ms (IQR 402.0 ms–418.6 ms) |
| jl_bpcells_rowmajor | 0.2 ms (IQR 0.2 ms–0.3 ms) | 0.8 ms (IQR 0.7 ms–0.8 ms) | 7.6 ms (IQR 7.2 ms–7.9 ms) | 70.9 ms (IQR 69.6 ms–71.8 ms) |
| raw_anndata_py | 38.2 ms (IQR 37.2 ms–39.3 ms) | 40.6 ms (IQR 37.9 ms–41.2 ms) | 51.1 ms (IQR 47.4 ms–51.6 ms) | 113.5 ms (IQR 110.6 ms–115.5 ms) |
| raw_h5py_py | 0.1 ms (IQR 0.1 ms–0.2 ms) | 1.2 ms (IQR 1.2 ms–1.2 ms) | 13.6 ms (IQR 12.8 ms–13.7 ms) | 134.5 ms (IQR 132.4 ms–135.3 ms) |
| py_scuba_sc | 24.7 ms (IQR 23.4 ms–29.1 ms) | 37.3 ms (IQR 33.2 ms–61.0 ms) | 214.4 ms (IQR 206.1 ms–237.8 ms) | 1.5 s (IQR 1.5 s–1.5 s) |
| r_seurat_bpcells [r] | 523.0 ms (IQR 519.0 ms–529.5 ms) | 595.0 ms (IQR 589.0 ms–616.0 ms) | 1.9 s (IQR 1.6 s–1.9 s) | 7.4 s (IQR 7.3 s–7.6 s) |
| r_seurat_bpcells_rowmajor [r] | 223.0 ms (IQR 219.0 ms–229.5 ms) | 288.0 ms (IQR 281.0 ms–296.0 ms) | 1.2 s (IQR 1.2 s–1.5 s) | 6.9 s (IQR 6.9 s–7.1 s) |
| r_sce_backed [r] | 50.0 ms (IQR 50.0 ms–54.0 ms) | 68.0 ms (IQR 67.0 ms–68.0 ms) | 264.0 ms (IQR 259.0 ms–273.0 ms) | 3.0 s (IQR 2.9 s–3.3 s) |
| raw_bpcells_r [r] | 328.0 ms (IQR 317.5 ms–330.5 ms) | 328.0 ms (IQR 324.5 ms–333.0 ms) | 414.0 ms (IQR 409.5 ms–425.5 ms) | 781.0 ms (IQR 771.0 ms–1.1 s) |
| raw_bpcells_r_rowmajor [r] | 10.0 ms (IQR 10.0 ms–10.0 ms) | 13.0 ms (IQR 13.0 ms–14.0 ms) | 55.0 ms (IQR 52.0 ms–58.0 ms) | 726.0 ms (IQR 439.0 ms–740.5 ms) |
| r_seurat_bpcells [julia] | 316.0 ms (IQR 313.5 ms–319.0 ms) | 322.0 ms (IQR 319.0 ms–324.5 ms) | 432.0 ms (IQR 428.0 ms–438.5 ms) | 829.0 ms (IQR 773.0 ms–1.1 s) |
| r_seurat_bpcells_rowmajor [julia] | 6.0 ms (IQR 5.0 ms–6.0 ms) | 9.0 ms (IQR 9.0 ms–12.5 ms) | 82.0 ms (IQR 44.5 ms–83.5 ms) | 480.0 ms (IQR 344.5 ms–768.5 ms) |
| r_sce_backed [julia] | 19.0 ms (IQR 18.5 ms–25.0 ms) | 58.0 ms (IQR 57.0 ms–63.5 ms) | 485.0 ms (IQR 459.0 ms–501.0 ms) | 4.4 s (IQR 4.4 s–4.5 s) |

## 100k cells: metadata, reductions, cells (median)

| implementation | fetch_cells:batch:b3 | fetch_cells:cell_type:T_cell | metadata:cell_type | metadata:full | metadata:nCount_RNA | reduction:X_pca|1-2 | reduction:X_pca|1-50 | reduction:X_umap|1-2 | unique_values:cell_type | unique_values:sample |
|---|---|---|---|---|---|---|---|---|---|---|
| jl_h5ad_csc | 13.1 ms | 13.9 ms | 5.8 ms | 29.1 ms | 0.2 ms | 19.7 ms | 496.8 ms | 18.3 ms | 1.3 ms | 1.3 ms |
| jl_h5ad_csr | 12.6 ms | 13.4 ms | 5.7 ms | 29.8 ms | 0.2 ms | 18.9 ms | 512.7 ms | 17.5 ms | 1.3 ms | 1.3 ms |
| jl_bpcells | 11.3 ms | 12.2 ms | 4.2 ms | 1.4 ms | 0.1 ms | 0.2 ms | 3.6 ms | 0.2 ms | 0.0 ms | 0.0 ms |
| jl_bpcells_rowmajor | 12.8 ms | 13.8 ms | 4.5 ms | 1.5 ms | 0.1 ms | 0.2 ms | 3.8 ms | 0.2 ms | 0.0 ms | 0.0 ms |
| raw_anndata_py | 98.6 ms | 103.8 ms | 200.5 ms | 0.1 ms | 0.1 ms | 0.7 ms | 1.8 ms | 0.2 ms | 192.1 ms | 185.0 ms |
| py_scuba_sc | 24.2 ms | 35.3 ms | 225.9 ms | 88.5 ms | 24.3 ms | 45.8 ms | 599.4 ms | 44.9 ms | 1.6 ms | 1.5 ms |
| r_seurat_bpcells [r] | 7.0 ms | 12.0 ms | 16.0 ms | 0.0 ms | 3.0 ms | 6.0 ms |  | 5.0 ms | 3.0 ms | 3.0 ms |
| r_seurat_bpcells_rowmajor [r] | 7.0 ms | 11.0 ms | 15.0 ms | 1.0 ms | 3.0 ms | 5.0 ms |  | 4.0 ms | 4.0 ms | 3.0 ms |
| r_sce_backed [r] | 5.0 ms | 8.0 ms | 16.0 ms | 8.0 ms | 4.0 ms | 25.0 ms |  | 25.0 ms | 5.0 ms | 5.0 ms |
| r_seurat_bpcells [julia] | 19.0 ms | 21.0 ms | 26.0 ms | 188.0 ms | 8.0 ms | 16.0 ms | 198.0 ms | 16.0 ms | 25.0 ms | 25.0 ms |
| r_seurat_bpcells_rowmajor [julia] | 20.0 ms | 21.0 ms | 27.0 ms | 154.0 ms | 8.0 ms | 15.0 ms | 184.0 ms | 16.0 ms | 26.0 ms | 25.0 ms |
| r_sce_backed [julia] | 20.0 ms | 23.0 ms | 27.0 ms | 149.0 ms | 8.0 ms | 21.0 ms | 215.0 ms | 22.0 ms | 26.0 ms | 27.0 ms |

## 1M cells: feature fetch over all cells, warm (median, IQR)

| implementation | genes=1 | genes=10 | genes=100 | genes=1000 |
|---|---|---|---|---|
| jl_h5ad_csc | 4.6 ms (IQR 2.8 ms–8.2 ms) | 59.0 ms (IQR 51.9 ms–74.3 ms) | 961.8 ms (IQR 870.5 ms–1.1 s) | 8.7 s (IQR 8.7 s–9.1 s) |
| jl_h5ad_csr | 92.0 s (IQR 91.3 s–94.7 s) | NA | NA | NA |
| jl_bpcells | 3.2 s (IQR 3.2 s–3.2 s) | 3.3 s (IQR 3.3 s–3.3 s) | 3.8 s (IQR 3.8 s–3.8 s) | 5.0 s (IQR 4.9 s–5.0 s) |
| jl_bpcells_rowmajor | 1.9 ms (IQR 1.9 ms–2.0 ms) | 8.2 ms (IQR 7.7 ms–8.3 ms) | 86.7 ms (IQR 86.0 ms–87.3 ms) | 1.7 s (IQR 1.6 s–1.9 s) |
| raw_anndata_py | 497.4 ms (IQR 479.4 ms–500.7 ms) | 464.3 ms (IQR 459.3 ms–479.7 ms) | 504.2 ms (IQR 502.0 ms–519.1 ms) | 809.3 ms (IQR 807.3 ms–815.9 ms) |
| raw_h5py_py | 0.1 ms (IQR 0.1 ms–0.1 ms) | 0.9 ms (IQR 0.8 ms–1.3 ms) | 13.1 ms (IQR 10.4 ms–14.6 ms) | 127.5 ms (IQR 126.8 ms–127.8 ms) |
| py_scuba_sc | 157.9 ms (IQR 72.3 ms–181.0 ms) | 281.5 ms (IQR 254.9 ms–294.3 ms) | 1.5 s (IQR 1.4 s–1.6 s) | 12.9 s (IQR 12.1 s–13.2 s) |
| r_seurat_bpcells [r] | 5.5 s (IQR 5.5 s–5.6 s) | 6.5 s (IQR 6.2 s–6.6 s) | 14.2 s (IQR 14.0 s–14.3 s) | 70.6 s (IQR 69.9 s–73.0 s) |
| r_seurat_bpcells_rowmajor [r] | 2.3 s (IQR 2.3 s–2.4 s) | 3.4 s (IQR 3.3 s–3.4 s) | 10.7 s (IQR 10.6 s–10.9 s) | 67.5 s (IQR 67.1 s–67.9 s) |
| r_sce_backed [r] | 254.0 ms (IQR 251.0 ms–256.0 ms) | 438.0 ms (IQR 431.0 ms–450.0 ms) | 3.9 s (IQR 3.8 s–4.1 s) | 39.6 s (IQR 39.4 s–40.1 s) |
| raw_bpcells_r [r] | 3.1 s (IQR 3.0 s–3.1 s) | 3.1 s (IQR 3.1 s–3.2 s) | 4.0 s (IQR 4.0 s–4.0 s) | 8.2 s (IQR 7.8 s–8.2 s) |
| raw_bpcells_r_rowmajor [r] | 85.0 ms (IQR 82.0 ms–87.5 ms) | 121.0 ms (IQR 115.5 ms–122.0 ms) | 479.0 ms (IQR 476.5 ms–490.0 ms) | 4.7 s (IQR 4.7 s–4.8 s) |
| r_seurat_bpcells [julia] | 3.3 s (IQR 3.2 s–3.3 s) | 3.4 s (IQR 3.4 s–3.4 s) | 4.1 s (IQR 4.1 s–4.3 s) | 8.0 s (IQR 7.8 s–8.1 s) |
| r_seurat_bpcells_rowmajor [julia] | 110.0 ms (IQR 107.0 ms–128.5 ms) | 202.0 ms (IQR 146.5 ms–212.0 ms) | 535.0 ms (IQR 436.0 ms–611.5 ms) | 4.7 s (IQR 4.6 s–4.9 s) |
| r_sce_backed [julia] | 181.0 ms (IQR 178.0 ms–187.0 ms) | 637.0 ms (IQR 628.5 ms–643.0 ms) | 4.5 s (IQR 4.1 s–4.6 s) | 40.6 s (IQR 40.2 s–40.9 s) |

## 1M cells: metadata, reductions, cells (median)

| implementation | fetch_cells:batch:b3 | fetch_cells:cell_type:T_cell | metadata:cell_type | metadata:full | metadata:nCount_RNA | reduction:X_pca|1-2 | reduction:X_pca|1-50 | reduction:X_umap|1-2 | unique_values:cell_type | unique_values:sample |
|---|---|---|---|---|---|---|---|---|---|---|
| jl_h5ad_csc | 147.0 ms | 157.7 ms | 56.5 ms | 389.6 ms | 1.7 ms | 204.9 ms | 4.9 s | 188.1 ms | 12.4 ms | 11.7 ms |
| jl_h5ad_csr | 141.8 ms | 156.0 ms | 57.5 ms | 446.1 ms | 1.7 ms | 203.1 ms | 5.0 s | 187.9 ms | 13.4 ms | 12.0 ms |
| jl_bpcells | 136.0 ms | 149.3 ms | 44.5 ms | 16.1 ms | 1.2 ms | 2.4 ms | 33.6 ms | 2.4 ms | 0.1 ms | 0.1 ms |
| jl_bpcells_rowmajor | 132.8 ms | 142.0 ms | 44.8 ms | 15.6 ms | 1.1 ms | 2.3 ms | 35.2 ms | 2.3 ms | 0.1 ms | 0.1 ms |
| raw_anndata_py | 55.3 ms | 85.2 ms | 130.3 ms | 0.3 ms | 0.3 ms | 12.4 ms | 44.9 ms | 0.6 ms | 108.9 ms | 78.1 ms |
| py_scuba_sc | 236.1 ms | 331.5 ms | 212.4 ms | 954.4 ms | 66.6 ms | 397.9 ms | 5.5 s | 379.4 ms | 14.2 ms | 13.3 ms |
| r_seurat_bpcells [r] | 85.0 ms | 137.0 ms | 144.0 ms | 2.0 ms | 32.0 ms | 60.0 ms |  | 62.0 ms | 34.0 ms | 34.0 ms |
| r_seurat_bpcells_rowmajor [r] | 88.0 ms | 133.0 ms | 140.0 ms | 2.0 ms | 28.0 ms | 56.0 ms |  | 53.0 ms | 33.0 ms | 35.0 ms |
| r_sce_backed [r] | 68.0 ms | 99.0 ms | 140.0 ms | 44.0 ms | 30.0 ms | 123.0 ms |  | 133.0 ms | 33.0 ms | 34.0 ms |
| r_seurat_bpcells [julia] | 236.0 ms | 291.0 ms | 340.0 ms | 1.8 s | 133.0 ms | 243.0 ms | 1.8 s | 232.0 ms | 401.0 ms | 429.0 ms |
| r_seurat_bpcells_rowmajor [julia] | 233.0 ms | 279.0 ms | 336.0 ms | 1.7 s | 135.0 ms | 252.0 ms | 1.7 s | 301.0 ms | 401.0 ms | 415.0 ms |
| r_sce_backed [julia] | 230.0 ms | 291.0 ms | 337.0 ms | 1.8 s | 139.0 ms | 259.0 ms | 1.9 s | 257.0 ms | 412.0 ms | 429.0 ms |

## Startup, open and first-call latency (median)

| implementation | scale | first_call | open | startup |
|---|---|---|---|---|
| jl_h5ad_csc | 100k | 167.7 ms | 400.4 ms | 3.9 s |
| jl_h5ad_csc | 10k | 146.5 ms | 362.5 ms | 1.8 s |
| jl_h5ad_csc | 1M | 336.3 ms | 714.0 ms | 1.6 s |
| jl_h5ad_csr | 100k | 30.6 ms | 30.0 ms |  |
| jl_h5ad_csr | 10k | 5.2 ms | 8.8 ms |  |
| jl_h5ad_csr | 1M | 326.5 ms | 259.7 ms |  |
| jl_bpcells | 100k | 12.3 ms | 2.0 s |  |
| jl_bpcells | 10k | 31.2 ms | 1.3 s | 1.2 s |
| jl_bpcells | 1M | 151.1 ms | 7.9 s |  |
| jl_bpcells_rowmajor | 100k | 29.2 ms | 2.0 s | 2.0 s |
| jl_bpcells_rowmajor | 10k | 29.1 ms | 1.5 s | 1.2 s |
| jl_bpcells_rowmajor | 1M | 29.0 ms | 7.8 s | 1.2 s |
| r_seurat_bpcells [r] | 100k | 280.0 ms | 712.0 ms | 4.8 s |
| r_seurat_bpcells [r] | 10k | 29.0 ms | 140.0 ms | 5.6 s |
| r_seurat_bpcells [r] | 1M | 2.6 s | 5.3 s | 5.7 s |
| r_seurat_bpcells_rowmajor [r] | 100k | 19.0 ms | 780.0 ms | 4.7 s |
| r_seurat_bpcells_rowmajor [r] | 10k | 5.0 ms | 152.0 ms | 4.8 s |
| r_seurat_bpcells_rowmajor [r] | 1M | 199.0 ms | 5.3 s | 4.8 s |
| r_sce_backed [r] | 100k | 35.0 ms | 215.0 ms |  |
| r_sce_backed [r] | 10k | 28.0 ms | 32.0 ms |  |
| r_sce_backed [r] | 1M | 134.0 ms | 2.3 s |  |
| raw_bpcells_r [r] | 100k | 327.0 ms | 24.0 ms |  |
| raw_bpcells_r [r] | 10k | 37.0 ms | 9.0 ms |  |
| raw_bpcells_r [r] | 1M | 3.1 s | 213.0 ms |  |
| raw_bpcells_r_rowmajor [r] | 100k | 12.0 ms | 24.0 ms |  |
| raw_bpcells_r_rowmajor [r] | 10k | 4.0 ms | 11.0 ms |  |
| raw_bpcells_r_rowmajor [r] | 1M | 85.0 ms | 197.0 ms |  |
| r_seurat_bpcells [julia] | 100k | 314.5 ms | 748.0 ms | 6.7 s |
| r_seurat_bpcells [julia] | 10k | 91.0 ms | 149.0 ms | 4.6 s |
| r_seurat_bpcells [julia] | 1M | 1.4 s | 5.1 s | 4.7 s |
| r_seurat_bpcells_rowmajor [julia] | 100k | 194.0 ms | 706.0 ms | 4.1 s |
| r_seurat_bpcells_rowmajor [julia] | 10k | 86.0 ms | 162.0 ms | 4.5 s |
| r_seurat_bpcells_rowmajor [julia] | 1M | 466.0 ms | 5.1 s | 4.6 s |
| r_sce_backed [julia] | 100k | 117.0 ms | 214.0 ms |  |
| r_sce_backed [julia] | 10k | 47.0 ms | 32.0 ms |  |
| r_sce_backed [julia] | 1M | 462.0 ms | 2.4 s |  |
| raw_anndata_py | 100k | 51.1 ms | 780.0 ms | 1.3 s |
| raw_anndata_py | 10k | 9.4 ms | 545.5 ms | 1.3 s |
| raw_anndata_py | 1M | 96.9 ms | 3.2 s | 1.6 s |
| raw_h5py_py | 100k | 0.3 ms | 49.9 ms |  |
| raw_h5py_py | 10k | 0.2 ms | 14.7 ms |  |
| raw_h5py_py | 1M | 0.2 ms | 403.2 ms |  |
| py_scuba_sc | 100k | 394.8 ms | 5.8 s |  |
| py_scuba_sc | 10k | 367.0 ms | 5.4 s | 35.1 ms |
| py_scuba_sc | 1M | 573.4 ms | 6.7 s |  |

## Skipped and failed cells

| implementation | engine | scale | status | note | rows |
|---|---|---|---|---|---|
| jl_h5ad_csr | native | 1M | skipped | rule: impl=jl_h5ad_csr scale=1M n_features>=10 | 63 |
| jl_bpcells | native | 100k | error | BoundsError: attempt to access 1-element Vector{SubString{String}} at index [2] |  3 |
| jl_bpcells_rowmajor | native | 100k | error | BoundsError: attempt to access 1-element Vector{SubString{String}} at index [2] |  3 |
| jl_h5ad_csc | native | 100k | error | BoundsError: attempt to access 1-element Vector{SubString{String}} at index [2] |  3 |
| r_sce_backed | r | 100k | error | NA/NaN argument |  3 |
| r_seurat_bpcells | julia | 100k | error | NA/NaN argument |  3 |
| r_seurat_bpcells | r | 100k | error | NA/NaN argument |  3 |
| r_seurat_bpcells_rowmajor | julia | 100k | error | NA/NaN argument |  3 |
| r_seurat_bpcells_rowmajor | r | 100k | error | NA/NaN argument |  3 |
| raw_anndata_py | none | 100k | error | not enough values to unpack (expected 2, got 1) |  3 |
| jl_bpcells | native | 1M | error | BoundsError: attempt to access 1-element Vector{SubString{String}} at index [2] |  2 |
| jl_bpcells_rowmajor | native | 1M | error | BoundsError: attempt to access 1-element Vector{SubString{String}} at index [2] |  2 |
| jl_h5ad_csc | native | 1M | error | BoundsError: attempt to access 1-element Vector{SubString{String}} at index [2] |  2 |
| r_sce_backed | r | 1M | error | NA/NaN argument |  2 |
| r_seurat_bpcells | julia | 1M | error | NA/NaN argument |  2 |
| r_seurat_bpcells | r | 1M | error | NA/NaN argument |  2 |
| r_seurat_bpcells_rowmajor | julia | 1M | error | NA/NaN argument |  2 |
| r_seurat_bpcells_rowmajor | r | 1M | error | NA/NaN argument |  2 |
| raw_anndata_py | none | 1M | error | not enough values to unpack (expected 2, got 1) |  2 |
| r_anndata_backed | r | 100k | skipped | R anndata 0.8.0 returns python-class objects; install remotes::install_version('anndata', '0.7.5.6') |  1 |
| r_anndata_backed | r | 10k | skipped | R anndata 0.8.0 returns python-class objects; install remotes::install_version('anndata', '0.7.5.6') |  1 |
| r_anndata_backed | r | 1M | skipped | R anndata 0.8.0 returns python-class objects; install remotes::install_version('anndata', '0.7.5.6') |  1 |
| r_sce_backed | r | 100k | skipped | legacy fetch_reduction supports exactly two dims |  1 |
| r_sce_backed | r | 10k | skipped | legacy fetch_reduction supports exactly two dims |  1 |
| r_sce_backed | r | 1M | skipped | legacy fetch_reduction supports exactly two dims |  1 |
| r_seurat_bpcells | r | 100k | skipped | legacy fetch_reduction supports exactly two dims |  1 |
| r_seurat_bpcells | r | 10k | skipped | legacy fetch_reduction supports exactly two dims |  1 |
| r_seurat_bpcells | r | 1M | skipped | legacy fetch_reduction supports exactly two dims |  1 |
| r_seurat_bpcells_rowmajor | r | 100k | skipped | legacy fetch_reduction supports exactly two dims |  1 |
| r_seurat_bpcells_rowmajor | r | 10k | skipped | legacy fetch_reduction supports exactly two dims |  1 |
| r_seurat_bpcells_rowmajor | r | 1M | skipped | legacy fetch_reduction supports exactly two dims |  1 |

## Checksum mismatches

none

## Figures

- `figures/cold_vs_warm.png`
- `figures/concurrency_p95.png`
- `figures/concurrency_throughput.png`
- `figures/features_100k.png`
- `figures/features_10k.png`
- `figures/features_1M.png`
- `figures/memory.png`
- `figures/other_100k.png`
- `figures/other_10k.png`
- `figures/other_1M.png`
- `figures/scaling_features_all.png`
- `figures/scaling_features_rand10.png`
