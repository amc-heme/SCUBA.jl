# Peak resident memory of the R process for a wide fetch (1000 genes x all cells) under the legacy
# engine and the Julia engine with double / single element type. Run under /usr/bin/time -l (macOS)
# or /usr/bin/time -v (Linux) and read "maximum resident set size":
#   Rscript bench/r/eltype_rss.R <scale> <r|julia> <double|single>
args <- commandArgs(trailingOnly = TRUE)
scale <- args[1]; engine <- args[2]; prec <- if (length(args) >= 3) args[3] else "double"
cfg <- RcppTOML::parseTOML(file.path(dirname(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))), "config", "bench.toml"))
suppressPackageStartupMessages(pkgload::load_all(path.expand(cfg$suite$r_scuba_path), quiet = TRUE))
options(SCUBA.engine = engine, SCUBA.julia_eltype = prec)
data_dir <- file.path(path.expand(Sys.getenv("SCUBA_BENCH_DATA", "~/scuba-bench-data")), scale)
obj <- readRDS(file.path(data_dir, "seurat_bpcells_rowmajor.rds"))
sel <- jsonlite::fromJSON(file.path(data_dir, "selection.json"))
genes <- paste0("rna_", sel$gene_sets[["1000:0"]])
if (engine == "julia") scuba_julia_available()
invisible(fetch_data(obj, genes[1]))                       # warm
t <- system.time(d <- fetch_data(obj, genes))
rss <- as.numeric(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE)) / 1024
cat(sprintf("scale=%s engine=%s eltype=%s cells=%d genes=%d wall=%.1fs rss_after_mb=%.0f\n", scale, engine, prec, nrow(d), ncol(d), t[["elapsed"]], rss))
