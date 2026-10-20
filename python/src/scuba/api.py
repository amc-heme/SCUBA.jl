"""Public API: pandas-flavoured wrappers around SCUBA.jl."""
from __future__ import annotations

import os
from typing import Any, Sequence

import numpy as np
import pandas as pd

from . import _jl


class Dataset:
    """Handle on a Julia-side SCUBA source (an opened file, or a wrapped live object)."""

    def __init__(self, jl_source, label: str = ""):
        self._src = jl_source
        self.label = label

    def __repr__(self) -> str:
        return f"scuba.Dataset({self.label or str(self._src)})"

    # convenience methods mirroring the module-level API
    def fetch_data(self, vars, layer=None, cells=None, strict=False):
        return fetch_data(self, vars, layer=layer, cells=cells, strict=strict)

    def fetch_metadata(self, vars=None, cells=None, full_table=False, return_class="dataframe"):
        return fetch_metadata(self, vars, cells=cells, full_table=full_table, return_class=return_class)

    def fetch_reduction(self, reduction=None, cells=None, dims=(1, 2)):
        return fetch_reduction(self, reduction, cells=cells, dims=dims)


def init(threads: int | None = None, sysimage: str | None = None) -> None:
    """Start the Julia runtime now (it otherwise starts on the first call)."""
    if threads is not None:
        os.environ.setdefault("PYTHON_JULIACALL_THREADS", str(threads))
    if sysimage is not None:
        os.environ.setdefault("PYTHON_JULIACALL_SYSIMAGE", sysimage)
    _jl.jl()


def open(path: str | os.PathLike) -> Dataset:
    """Open an on-disk container with SCUBA.jl's native readers (.h5ad, zarr store, BPCells directory)."""
    p = os.fspath(path)
    S = _jl.S()
    if os.path.isdir(p):
        if os.path.exists(os.path.join(p, "zarr.json")) or os.path.exists(os.path.join(p, ".zgroup")):
            return Dataset(_jl.call(S.open_zarr, p), p)
        if os.path.exists(os.path.join(p, "version")) and os.path.exists(os.path.join(p, "idxptr")):
            raise ValueError("a BPCells directory is a single matrix; build a dataset around it from Julia (SCUBA.bpcells_modality)")
        raise ValueError(f"unrecognised directory container: {p}")
    if p.endswith(".h5ad") or p.endswith(".h5"):
        return Dataset(_jl.call(S.open_h5ad, p), p)
    raise ValueError(f"unrecognised container: {p}")


def _wrap(obj) -> Dataset:
    if isinstance(obj, Dataset):
        return obj
    if isinstance(obj, (str, os.PathLike)):
        return open(obj)
    # live anndata.AnnData (duck-typed to avoid importing anndata)
    if hasattr(obj, "obs_names") and hasattr(obj, "var_names"):
        S = _jl.S()
        return Dataset(_jl.call(S.pysource, obj), "anndata.AnnData")
    raise TypeError(f"unsupported object of type {type(obj).__name__}")


def _vars(vars) -> list[str]:
    if isinstance(vars, str):
        return [vars]
    return [str(v) for v in vars]


def _cells(cells):
    if cells is None:
        return None
    return [str(c) for c in cells]


def _to_pandas(df) -> pd.DataFrame:
    """Julia DataFrame -> pandas, built on the Julia side (cell IDs as index, categoricals preserved)."""
    return _jl.glue().to_pandas(df)


def _to_list(v) -> list:
    return [x.item() if hasattr(x, "item") else (str(x) if not isinstance(x, (int, float, bool)) else x) for x in v]


def fetch_data(obj, vars: str | Sequence[str], layer: str | None = None, cells=None, strict: bool = False) -> pd.DataFrame:
    """Fetch features, metadata columns and reduction coordinates into one DataFrame (index = cell IDs)."""
    ds = _wrap(obj)
    S = _jl.S()
    df = _jl.call(S.fetch_data, ds._src, _vars(vars), layer=layer, cells=_cells(cells), strict=strict)
    return _to_pandas(df)


def fetch_feature(obj, features, modality: str | None = None, layer: str | None = None, cells=None) -> pd.DataFrame:
    ds = _wrap(obj)
    S = _jl.S()
    df = _jl.call(S.fetch_feature, ds._src, _vars(features), modality=modality, layer=layer, cells=_cells(cells))
    return _to_pandas(df)


def fetch_metadata(obj, vars=None, cells=None, full_table: bool = False, return_class: str = "dataframe"):
    """Metadata columns. `vars=None` or `full_table=True` returns the whole table; `return_class="vector"` a Series."""
    ds = _wrap(obj)
    S = _jl.S()
    if return_class not in ("dataframe", "vector"):
        raise ValueError("return_class must be 'dataframe' or 'vector'")
    if vars is None or full_table:
        return _to_pandas(_jl.call(S.fetch_metadata, ds._src, cells=_cells(cells)))
    df = _to_pandas(_jl.call(S.fetch_metadata, ds._src, _vars(vars), cells=_cells(cells)))
    if return_class == "vector":
        if df.shape[1] != 1:
            raise ValueError("return_class='vector' requires exactly one var")
        return df.iloc[:, 0]
    return df


def fetch_reduction(obj, reduction: str | None = None, cells=None, dims=(1, 2)) -> pd.DataFrame:
    ds = _wrap(obj)
    S = _jl.S()
    if reduction is None:
        df = _jl.call(S.fetch_reduction, ds._src, cells=_cells(cells), dims=[int(d) for d in dims])
    else:
        df = _jl.call(S.fetch_reduction, ds._src, reduction, cells=_cells(cells), dims=[int(d) for d in dims])
    return _to_pandas(df)


def fetch_cells(obj, meta_var: str, meta_levels) -> list[str]:
    ds = _wrap(obj)
    S = _jl.S()
    levels = [meta_levels] if isinstance(meta_levels, (str, int, float, bool)) else list(meta_levels)
    return [str(c) for c in _jl.call(S.fetch_cells, ds._src, meta_var, [str(l) for l in levels])]


def get_all_cells(obj) -> list[str]:
    return [str(c) for c in _jl.call(_jl.S().all_cells, _wrap(obj)._src)]


def meta_varnames(obj) -> list[str]:
    return [str(c) for c in _jl.call(_jl.S().meta_varnames, _wrap(obj)._src)]


def all_keys(obj) -> dict[str, str]:
    d = _jl.call(_jl.S().all_keys, _wrap(obj)._src)
    return {str(k): str(v) for k, v in d.items()}


def features_in_assay(obj, assay: str | None = None) -> list[str]:
    S = _jl.S()
    src = _wrap(obj)._src
    res = _jl.call(S.features_in_modality, src) if assay is None else _jl.call(S.features_in_modality, src, assay)
    return [str(f) for f in res]


def reduction_dimnames(obj, reduction: str, dims) -> list[str]:
    return [str(n) for n in _jl.call(_jl.S().reduction_dimnames, _wrap(obj)._src, reduction, [int(d) for d in dims])]


def default_layer(obj) -> str:
    return str(_jl.call(_jl.S().default_layer, _wrap(obj)._src))


def default_reduction(obj) -> str:
    return str(_jl.call(_jl.S().default_reduction, _wrap(obj)._src))


def unique_values(obj, var: str) -> list:
    v = _jl.call(_jl.S().unique_values, _wrap(obj)._src, var)
    return [str(x) for x in v]


def n_unique(obj, var: str) -> int:
    return int(_jl.call(_jl.S().n_unique, _wrap(obj)._src, var))
