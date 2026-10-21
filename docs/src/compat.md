# Compatibility with SCUBA (R)

SCUBA.jl keeps the R package's API shape (`fetch_data`, `fetch_metadata`, `fetch_reduction`,
`fetch_cells`, `get_all_cells`/`all_cells`, `meta_varnames`, `all_keys`, `features_in_assay`/
`features_in_modality`, `reduction_dimnames`, `default_layer`, `default_reduction`,
`unique_values`, `n_unique`) and its messages, but implements the variable-resolution rules once.
The R package 1.2.2 behaved slightly differently per object class; SCUBA.jl normalises these
cases. The golden tests (`fixtures/cases.json`) flag each divergence.

| Behaviour | SCUBA R 1.2.2 | SCUBA.jl |
|---|---|---|
| Key matching | SCE matched `^key` without the underscore (`AB` matched `ABC_x`); `fetch_feature` on Seurat was case-insensitive | keys match as `key_`, case-sensitive; not-found errors suggest case variants |
| `X` vs `X_umap` (AnnData) | vars containing another key anywhere were excluded from `X` | longest matching key whose remainder exists wins; ambiguity warns |
| Duplicates | Seurat path silent and returned both spellings; SCE/AnnData warned | always collapsed by resolved target, with a warning |
| Output names (AnnData) | `obs_` stripped always, `X_` stripped only when the user omitted it | keyed names kept as written, `obs_` stripped |
| Single reduction dimension on SCE | column dropped or misnamed | `UMAP_1` works |
| Reduction column names on SCE | matrix column names when present | always `key_dim` |
| Unknown cells | class-dependent | `CellsNotFound` error |
| `fetch_feature` names | keyed names from any assay (case-insensitive), grouped per assay; a *bare* name that is also a metadata column resolved metadata-first through `fetch_data` on non-BPCells assays | keyed names from any modality (case-sensitive), grouped per modality, requested order and names kept; bare names are only searched in the main modality (the R shim keeps colliding bare names on the legacy method) |
| Large queries | warning at ≥ 1000 vars | warning when the dense result exceeds `SCUBA.MAX_DENSE_BYTES[]` (1 GiB) |
| Layer missing in a single-layer modality (e.g. `obsm` data frame) | ignored (AnnData) or error (SCE) | the only layer is used, with a `:layer_substituted` warning |
| `fetch_metadata(return_class = "vector")` | argument | method `fetch_metadata(src, var)` |
| `update_object_metadata` | returns a modified copy | `set_metadata!` mutates in place; R/Seurat bridges are read-only (assign in R) |

Everything else is identical: metadata beats a same-named feature (with a warning), an unkeyed
feature found in exactly one other modality is returned renamed `key_feature` with a warning,
partial misses warn and total misses error with the same wording, `fetch_metadata` errors on any
missing variable, `fetch_cells` rejects numeric variables, `default_reduction` prefers umap, then
tsne, then pca.
