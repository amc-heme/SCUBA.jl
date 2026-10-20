# Fetching data

## Keys

Every modality (Seurat assay, SCE experiment, AnnData `X` or `obsm` data frame) and every
reduction has a *key*. A variable is addressed as `key_name`:

| Source | modality keys | reduction keys | examples |
|---|---|---|---|
| AnnData (`.h5ad`, zarr, Python) | `X`, each `obsm` data frame (e.g. `protein`), `raw` | each numeric `obsm` array (`X_umap`, `X_pca`) | `X_GAPDH`, `protein_CD9-AB`, `X_umap_1` |
| Seurat (RCall) | assay keys without the underscore (`rna`, `ab`); assay names are aliases | reduction keys (`UMAP`, `PC`); names (`umap`) are aliases | `rna_GAPDH`, `ab_CD9-AB`, `UMAP_1` |
| SingleCellExperiment (RCall) | `mainExpName`, `altExpNames` | `reducedDimNames` | `RNA_GAPDH`, `AB_CD9-AB`, `UMAP_1` |
| `MemoryAnnData` | the modality keys you chose | the reduction keys you chose | |

Metadata columns are addressed by name, or as `obs_name`. `all_keys(src)` lists the keys.

## Resolution rules

`fetch_data(src, vars)` resolves each var with one rule set, identical for every backend:

1. **Keyed names first.** A key matches only as `key_` (with the underscore), exactly and
   case-sensitively. When several keys match (`X` and `X_umap` for `X_umap_1`), the longest key
   whose remainder exists wins; if more than one exists an `:ambiguous_var` warning is emitted.
2. **Unkeyed names** are looked up as metadata, then as a feature of the main modality, then as a
   feature of exactly one other modality. In the last case the column is renamed `key_name` and
   an `:alt_modality_renamed` warning is emitted. A name present in several other modalities is
   ambiguous and treated as missing.
3. **Metadata shadows features.** If an unkeyed name is both a metadata column and a feature of
   the main modality, the metadata is returned with a `:shadowed_feature` warning; use the keyed
   form for the feature.
4. **Duplicates** collapse by resolved target, so `GAPDH` and `RNA_GAPDH` yield one column.
5. **Reductions** take an integer dimension (`UMAP_2`); out-of-range dimensions are not found.
6. **Missing variables**: if none resolve, `VarsNotFound` is thrown; otherwise a
   `:some_vars_not_found` warning lists them (`strict = true` turns it into the error).
7. Columns come back in input order. Keyed names are kept as written, `obs_` is stripped.

## Result tables

Results are `DataFrame`s whose first column holds the cell identifiers. Its name is normally
`cell` (see [`cell_column`](@ref) when a fetched variable is itself called `cell`). Categorical
metadata is a `CategoricalVector` with level order preserved; missing values are `missing`.
Each column carries `colmetadata(df, col, "scuba")`, a named tuple `(kind, target, layer, input)`
recording where it came from.

## Other accessors

- [`fetch_feature`](@ref) — features from one modality only (fast path for plotting).
- [`fetch_matrix`](@ref) — a dense or sparse cells × features matrix.
- [`fetch_metadata`](@ref), [`fetch_cells`](@ref), [`unique_values`](@ref), [`n_unique`](@ref).
- [`fetch_reduction`](@ref), [`default_reduction`](@ref), [`reduction_dimnames`](@ref).
- [`snapshot`](@ref) — materialise any source in memory for repeated queries.
- [`with_collected_warnings`](@ref) — run a call and collect its structured warnings.
