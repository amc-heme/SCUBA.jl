import os
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
FIXTURES = pathlib.Path(os.environ.get("SCUBA_FIXTURES_DIR", REPO / "fixtures" / "data"))

# Use the development checkout of SCUBA.jl and an isolated Julia depot/env for the tests.
os.environ.setdefault("PYTHON_JULIAPKG_PROJECT", str(REPO / "python" / ".juliapkg-test"))
os.environ.setdefault("JULIA_CONDAPKG_BACKEND", "Null")


@pytest.fixture(scope="session")
def h5ad_path():
    p = FIXTURES / "AML_h5ad.h5ad"
    if not p.exists():
        pytest.skip(f"fixture not found: {p}")
    return str(p)


@pytest.fixture(scope="session")
def adata(h5ad_path):
    ad = pytest.importorskip("anndata")
    return ad.read_h5ad(h5ad_path)
