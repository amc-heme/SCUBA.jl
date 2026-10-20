# Result row schema

Every harness appends one CSV row per repetition (never pre-aggregated) to
`results/raw/<run_id>/<language>_<implementation>_<scale>_<cache_state>.csv`, flushing after each row.
`tools/merge_results.py` concatenates all raw files into `results/results.parquet` and `results/results.csv`.

| column | type | meaning |
|---|---|---|
| suite_version | str | `[suite].version` from bench.toml |
| run_id | str | run identifier (timestamp unless configured) |
| timestamp | str | ISO-8601 UTC when the row was written |
| host | str | `hostname` |
| language | str | `julia`, `r`, `python` |
| implementation | str | see list below |
| engine | str | `native` (Julia / Python), `r` or `julia` (R package engine) |
| storage | str | `h5ad_csc`, `h5ad_csr`, `bpcells`, `memory`, `sce_h5sparse`, `seurat_rds` |
| scale | str | `10k`, `100k`, `1M`, `3M` |
| workload | str | `feature`, `matrix_sparse`, `metadata`, `reduction`, `fetch_cells`, `unique_values`, `julia_rtt`, `open`, `startup`, `snapshot` |
| n_features | int | number of genes requested (0 when not applicable) |
| gene_set | int | gene set replicate index (0-based; -1 when not applicable) |
| cell_subset | str | `all`, `rand1`, `rand10`, `contig1`, `contig10` (empty when not applicable) |
| layer | str | `data`, `counts`, or empty |
| n_cells_out | int | rows in the returned table |
| cache_state | str | `warm`, `cold`, `cold_unverified` |
| phase | str | `startup`, `open`, `first_call`, `fetch`, `concurrent` |
| rep | int | repetition index (0 = discarded warm-up / first call) |
| threads | int | Julia threads or R/Python worker count |
| wall_ms | float | wall-clock time of the call |
| cpu_ms | float | user + system CPU time of the process during the call |
| open_ms | float | time to open the source (phase = open) or 0 |
| rss_mb | float | process RSS after the call |
| maxrss_mb | float | process high-water RSS so far |
| alloc_mb | float | memory allocated during the call (Julia bytes, R gc delta, Python tracemalloc) or NA |
| gc_ms | float | time in garbage collection during the call (Julia) or NA |
| result_checksum | str | value computed from the returned data (see expected.json) |
| checksum_ok | bool | matches expected.json within rtol 1e-5 |
| status | str | `ok`, `skipped`, `skipped_extrapolated`, `timeout`, `error`, `cold_unverified` |
| note | str | free text: tunables, error message, skip reason |

## Implementations

| implementation | language | engine | storage |
|---|---|---|---|
| jl_h5ad_csc | julia | native | h5ad_csc |
| jl_h5ad_csr | julia | native | h5ad_csr |
| jl_bpcells | julia | native | bpcells |
| jl_snapshot | julia | native | memory |
| r_seurat_bpcells | r | r / julia | bpcells (Seurat v5 object) |
| r_sce_backed | r | r / julia | sce_h5sparse (SingleCellExperiment over HDF5Array::H5SparseMatrix) |
| r_anndata_backed | r | r | h5ad_csc (reticulate, backed) |
| r_seurat_inmem | r | r | memory (dgCMatrix Seurat object, opt-in, <= 1M) |
| raw_bpcells_r | r | none | bpcells (`as.matrix(mat[genes, cells])`) |
| raw_seurat_fetchdata | r | none | memory (`SeuratObject::FetchData`, opt-in) |
| raw_anndata_py | python | none | h5ad_csc (`read_h5ad(backed="r")`) |
| raw_h5py_py | python | none | h5ad_csc (direct indptr slicing) |
| py_scuba_sc | python | native | h5ad_csc (scuba-sc through juliacall) |

## Checksums (expected.json)

- feature / matrix_sparse: `sum` of the fetched block in float64, `nnz`, and the first nonzero as `[cell_id, gene, value]`; the harness reports `sum` as `result_checksum` and compares with rtol 1e-5.
- metadata: `"<nrow>x<ncol>:<sum nCount_RNA>"` for the full table; `"<n>:<sum>"` for numeric columns; `"<n>:<sha1 of joined values>"` for categorical columns.
- reduction: sum of the requested coordinates in float64.
- fetch_cells: `"<count>:<sha1 of sorted ids joined by newline>"`.
- unique_values: `"<n>:<sha1 of sorted levels joined by newline>"`.
