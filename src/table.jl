# Assembly of the result DataFrame and the cell-column convention.

const CELL_COLUMN_KEY = "scuba_cell_column"
const COLUMN_META_KEY = "scuba"

"""
    cell_column(df) -> String

Name of the column holding cell identifiers in a table returned by SCUBA (`"cell"` unless
a fetched variable had that name).
"""
cell_column(df::AbstractDataFrame) = String(metadata(df, CELL_COLUMN_KEY, "cell"))

"""`cells(df)` — the cell identifier column of a SCUBA result table."""
cells(df::AbstractDataFrame) = df[!, cell_column(df)]

function _pick_cell_column_name(colnames)
    for cand in ("cell", "cell_id", "scuba_cell")
        cand in colnames || return cand
    end
    i = 1
    while string("cell_", i) in colnames
        i += 1
    end
    return string("cell_", i)
end

"""
    assemble_table(cell_names, columns, resolved) -> DataFrame

Build the result table: cell IDs first, then `columns` (name => vector) in order. Each
column carries `colmetadata(df, col, "scuba")` describing where the data came from.
"""
function assemble_table(cell_names::AbstractVector{<:AbstractString},
                        columns::AbstractVector{<:Pair{String}},
                        resolved::Union{Nothing,AbstractVector{ResolvedVar}} = nothing)
    name = _pick_cell_column_name(first.(columns))
    df = DataFrame(name => collect(String, cell_names); copycols = false)
    for (n, v) in columns
        df[!, n] = v
    end
    metadata!(df, CELL_COLUMN_KEY, name; style = :note)
    if resolved !== nothing
        for (k, rv) in enumerate(resolved)
            col = columns[k].first
            colmetadata!(df, col, COLUMN_META_KEY,
                         (kind = rv.kind, target = rv.target, layer = rv.layer, input = rv.input);
                         style = :note)
        end
    end
    return df
end

"""Table without the cell column (for callers that want plain data)."""
drop_cells(df::AbstractDataFrame) = select(df, Not(cell_column(df)))
