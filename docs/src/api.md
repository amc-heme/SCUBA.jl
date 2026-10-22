# API reference

## Accessors

```@docs
fetch_data
fetch_feature
fetch_matrix
eachchunk
fetch_metadata
fetch_reduction
fetch_cells
all_cells
meta_varnames
all_keys
features_in_modality
reduction_dimnames
default_reduction
unique_values
n_unique
set_metadata!
update_object_metadata!
snapshot
resolve_vars
cell_column
cells
drop_cells
with_collected_warnings
```

## Sources

```@docs
AbstractSingleCellSource
MemoryAnnData
Modality
SCUBA.LazyMatrix
H5ADSource
BPCellsMatrix
SCUBA.bpcells_modality
rsource
pysource
open_zarr
```

## Backend interface

```@docs
SCUBA.cell_ids
SCUBA.modality_keys
SCUBA.default_layer
SCUBA.expression
SCUBA.metadata_column
SCUBA.reduction
SCUBA.resolve_layer
SCUBA.canonical_modality
SCUBA.canonical_reduction
```

## Errors and warnings

```@docs
VarsNotFound
MetadataNotFound
SCUBAWarning
SCUBA.scuba_warn
```

## Fixtures

```@docs
fixtures_dir
fixture_path
```
