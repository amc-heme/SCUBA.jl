# A small synthetic MemoryAnnData with known values, used for unit tests of the resolution
# rules and the API. Layout mirrors the AML fixtures: main modality "RNA" (features g1..g5),
# alt modality "AB" (p1..p3, and "shared" which also exists in RNA), reductions UMAP (2 dims)
# and PCA (3 dims), metadata with a categorical, a numeric, a bool and a column named like a gene.
function synthetic_source()
    cells = ["c$i" for i in 1:6]
    rna_features = ["g1", "g2", "g3", "shared", "g5"]
    ab_features = ["p1", "p2", "shared", "obs_trap"]
    rna_data = sparse(Float64[i + 10j for i in 1:6, j in 1:5])
    rna_counts = sparse(Float64[i * j for i in 1:6, j in 1:5])
    ab_data = Float64[100 + i + j for i in 1:6, j in 1:4]
    obs = DataFrame(
        group = categorical(["a", "a", "b", "b", "c", "c"]; levels = ["c", "b", "a"]),
        score = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        flag = [true, false, true, false, true, false],
        g1 = ["m1", "m2", "m3", "m4", "m5", "m6"],          # metadata shadowing a RNA feature
        withmissing = [1, missing, 3, missing, 5, 6],
    )
    reductions = OrderedDict("UMAP" => Float64[i + 0.5j for i in 1:6, j in 1:2],
                             "PCA" => Float64[-i * j for i in 1:6, j in 1:3])
    MemoryAnnData(; cells, obs,
                  modalities = OrderedDict("RNA" => Modality(rna_features, OrderedDict("data" => rna_data, "counts" => rna_counts)),
                                           "AB" => Modality(ab_features, OrderedDict("data" => ab_data))),
                  reductions, default_layer = "data")
end
