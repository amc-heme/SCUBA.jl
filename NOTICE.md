# Notices

SCUBA.jl is a Julia reimplementation of the R package SCUBA
(<https://github.com/amc-heme/SCUBA>, MIT, Copyright (c) 2023 SCUBA Authors).

The accessor API (`fetch_data` and friends) is modelled on the design of
`SeuratObject::FetchData()` from the Seurat project (MIT, Copyright (c) 2021
Seurat authors). No Seurat source code is included in SCUBA.jl.

The bundled test fixtures are a 250-cell down-sample of the CITE-seq dataset of
Triana et al. (2021), *Nature Immunology*, distributed by the authors on
Figshare. See `fixtures/README.md` for provenance.
