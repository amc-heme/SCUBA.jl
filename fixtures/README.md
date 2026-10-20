# Fixtures

All fixtures derive from one 250-cell down-sample (10 `condensed_cell_type` groups × 25 cells)
of the CITE-seq healthy blood and bone marrow reference of Triana et al. 2021 (*Nature
Immunology*; data distributed by the authors on Figshare). The original down-sampling is
documented in `amc-heme/SCUBA_Manuscript/Demo_Object_Generation.Rmd`; the four base forms are
copied from the R package SCUBA 1.2.2 (`data/` and `inst/extdata/`).

| File | Form | Produced by |
|---|---|---|
| `AML_h5ad.h5ad` | AnnData, CSC X + `layers/counts`, obsm `X_pca`, `X_umap`, `protein` (dataframe) | SCUBA R package |
| `AML_sce/` | SingleCellExperiment saved with `HDF5Array::saveHDF5SummarizedExperiment` | SCUBA R package |
| `AML_BPCells/{RNA,AB}_{counts,data}/` | BPCells matrix directories | SCUBA R package |
| `AML_Seurat.rda`, `AML_BPCells.rda` | Seurat objects (in-memory / BPCells layers) | SCUBA R package |
| `AML_h5ad_csr.h5ad`, `AML_h5ad_dense.h5ad`, `AML.zarr`, `AML.h5mu`, `AML_h5ad_edge.h5ad` | derived AnnData/MuData forms | `make_fixtures.py` |
| `build/goldens/` | expected outputs of SCUBA R 1.2.2 for `cases.json` | `make_goldens.R` |

Regenerate with `make -C fixtures` (needs R with SCUBA 1.2.2, and the `scuba-fixtures` conda env).
The working copies live in `fixtures/data/` (git-ignored) and are shipped to users as the lazy
`scuba_fixtures` artifact (`pack.jl`). Tests find them through `SCUBA.fixtures_dir()`.
