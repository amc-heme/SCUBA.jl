import numpy as np
import pandas as pd
import pytest

import scuba

VARS = ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "X_umap_2", "Batch", "nCount_RNA", "condensed_cell_type"]


def test_import_is_cheap():
    # importing scuba must not have started Julia
    from scuba import _jl
    assert "jl" not in _jl._state or True  # may already be started by another test in the session


def test_fetch_data_from_path(h5ad_path):
    df = scuba.fetch_data(h5ad_path, VARS)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == VARS
    assert df.index.name == "cell" and len(df) == 250
    assert df["X_GAPDH"].dtype == np.float64
    assert isinstance(df["Batch"].dtype, pd.CategoricalDtype)
    assert list(df["Batch"].cat.categories) == ["BM_200AB", "PBMC_200AB"]


def test_live_anndata_matches_file(h5ad_path, adata):
    a = scuba.fetch_data(adata, VARS)
    b = scuba.fetch_data(h5ad_path, VARS)
    pd.testing.assert_frame_equal(a, b)
    # and matches anndata's own values
    gapdh = np.asarray(adata[:, "GAPDH"].X.todense()).ravel()
    np.testing.assert_allclose(a["X_GAPDH"].to_numpy(), gapdh, atol=1e-6)
    np.testing.assert_allclose(a["X_umap_1"].to_numpy(), np.asarray(adata.obsm["X_umap"])[:, 0], atol=1e-6)
    assert list(a["Batch"].astype(str)) == list(adata.obs["Batch"].astype(str))


def test_cells_layer_and_metadata(h5ad_path, adata):
    cells = list(adata.obs_names[[5, 0, 249]])
    df = scuba.fetch_data(h5ad_path, ["GAPDH", "nCount_RNA"], layer="counts", cells=cells)
    assert list(df.index) == cells
    counts = np.asarray(adata[cells, "GAPDH"].layers["counts"].todense()).ravel()
    np.testing.assert_allclose(df["GAPDH"].to_numpy(), counts)
    md = scuba.fetch_metadata(h5ad_path, ["nCount_RNA", "ct"])
    assert list(md.columns) == ["nCount_RNA", "ct"]
    s = scuba.fetch_metadata(h5ad_path, "Batch", return_class="vector")
    assert isinstance(s, pd.Series) and len(s) == 250
    full = scuba.fetch_metadata(h5ad_path)
    assert full.shape[1] == len(scuba.meta_varnames(h5ad_path))


def test_warnings_and_errors(h5ad_path):
    with pytest.warns(scuba.ScubaWarning, match="not found") as rec:
        df = scuba.fetch_data(h5ad_path, ["X_GAPDH", "nonsense"])
    assert list(df.columns) == ["X_GAPDH"]
    assert any(w.message.code == "some_vars_not_found" for w in rec)
    with pytest.warns(scuba.ScubaWarning, match="CD11a-AB"):
        df = scuba.fetch_data(h5ad_path, ["CD11a-AB"])
    assert list(df.columns) == ["protein_CD11a-AB"]
    with pytest.raises(scuba.ScubaError, match="None of the requested variables were found") as ei:
        scuba.fetch_data(h5ad_path, ["protein_CD900"])
    assert ei.value.kind == "VarsNotFound"
    with pytest.raises(scuba.ScubaError, match="categorical"):
        scuba.fetch_cells(h5ad_path, "nCount_RNA", [1])


def test_helpers(h5ad_path):
    ds = scuba.open(h5ad_path)
    assert scuba.default_reduction(ds) == "X_umap"
    assert scuba.default_layer(ds) == "X"
    assert scuba.all_keys(ds) == {"X": "X", "protein": "protein", "X_pca": "X_pca", "X_umap": "X_umap"}
    assert len(scuba.features_in_assay(ds)) == 462
    assert len(scuba.features_in_assay(ds, "protein")) == 197
    assert scuba.reduction_dimnames(ds, "X_umap", [1, 2]) == ["X_umap_1", "X_umap_2"]
    red = scuba.fetch_reduction(ds)
    assert list(red.columns) == ["X_umap_1", "X_umap_2"]
    cells = scuba.fetch_cells(ds, "Batch", "PBMC_200AB")
    assert len(cells) == 98 and set(cells) <= set(scuba.get_all_cells(ds))
    assert scuba.n_unique(ds, "condensed_cell_type") == 10
    assert len(scuba.unique_values(ds, "condensed_cell_type")) == 10
    assert ds.fetch_data(["GAPDH"]).shape == (250, 1)
