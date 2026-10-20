"""Lazy Julia bootstrap. The runtime starts on first use; `scuba.init()` starts it eagerly."""
from __future__ import annotations

import os
import threading
import warnings

from ._errors import ScubaError, ScubaWarning

_lock = threading.Lock()
_state: dict = {}

# Julia-side glue: run a thunk, collect SCUBA warnings, return (result, warnings) as Python-friendly values.
_GLUE = r"""
module ScubaPy
using SCUBA, PythonCall, DataFrames, CategoricalArrays
const S = SCUBA

# Python lists arrive as PyList; SCUBA wants plain Julia vectors of strings / ints
tojl(x) = x
tojl(x::PyList) = isempty(x) ? String[] : (all(e -> e isa Integer, x) ? Int[pyconvert(Int, e) for e in x] : String[pyconvert(String, e) for e in x])

function call(f, args...; kwargs...)
    jargs = map(tojl, args)
    jkw = (; (k => tojl(v) for (k, v) in pairs(kwargs))...)
    result, ws = SCUBA.with_collected_warnings(() -> f(jargs...; jkw...))
    return result, [(string(w.code), w.message, w.vars) for w in ws]
end

to_pandas(df::DataFrame) = SCUBA.to_pandas(df)

errkind(e) = string(nameof(typeof(e)))

# Vector of strings -> Python list of str in one call (avoids per-element conversion from Python)
pystrings(v) = pylist([ismissing(x) ? pybuiltins.None : pystr(string(x)) for x in v])
end
"""


def _start():
    with _lock:
        if "jl" in _state:
            return _state["jl"]
        from juliacall import Main as jl  # noqa: PLC0415 - deliberate lazy import

        jl.seval("using SCUBA, PythonCall, DataFrames, CategoricalArrays")
        try:
            jl.seval("using Zarr")
        except Exception:  # Zarr is optional
            pass
        jl.seval(_GLUE)
        _state["jl"] = jl
        _state["S"] = jl.SCUBA
        _state["glue"] = jl.ScubaPy
        return jl


def jl():
    return _state.get("jl") or _start()


def S():
    jl()
    return _state["S"]


def glue():
    jl()
    return _state["glue"]


def call(fn, *args, **kwargs):
    """Call a Julia SCUBA function, re-emitting SCUBA warnings and mapping errors."""
    from juliacall import JuliaError  # noqa: PLC0415

    g = glue()
    try:
        result, ws = g.call(fn, *args, **kwargs)
    except JuliaError as e:  # pragma: no cover - exercised in tests through specific cases
        msg = str(e)
        kind = None
        try:
            kind = str(g.errkind(e.exception))
        except Exception:
            pass
        # strip Julia's "ErrorType: " style prefix if present
        if ": " in msg and msg.split(": ", 1)[0].isidentifier():
            msg = msg.split(": ", 1)[1]
        err = ScubaError(msg, kind)
        if os.environ.get("SCUBA_DEBUG"):
            raise err from e
        raise err from None
    for code, message, vars_ in ws:
        warnings.warn(ScubaWarning(str(message), str(code), [str(v) for v in vars_]), stacklevel=3)
    return result
