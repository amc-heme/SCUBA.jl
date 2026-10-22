"""scuba: Single Cell Unified Back-end API for Python, powered by SCUBA.jl.

Nothing Julia-related happens at import time; the Julia runtime starts on the first call.
"""
from ._errors import ScubaError, ScubaWarning
from .api import (
    Dataset,
    all_keys,
    default_layer,
    default_reduction,
    features_in_assay,
    fetch_cells,
    fetch_data,
    fetch_feature,
    fetch_matrix,
    fetch_metadata,
    fetch_reduction,
    get_all_cells,
    init,
    iter_chunks,
    meta_varnames,
    n_unique,
    open,
    reduction_dimnames,
    unique_values,
)

__all__ = [
    "Dataset", "ScubaError", "ScubaWarning", "all_keys", "default_layer", "default_reduction",
    "features_in_assay", "fetch_cells", "fetch_data", "fetch_feature", "fetch_matrix", "fetch_metadata",
    "fetch_reduction", "get_all_cells", "init", "iter_chunks", "meta_varnames", "n_unique", "open",
    "reduction_dimnames", "unique_values",
]
__version__ = "0.1.0"
