#!/usr/bin/env Rscript
# R benchmark harness: SCUBA R package (legacy engine or Julia engine) and raw R baselines.
#
#   Rscript bench/r/run.R --scale 10k --engine r|julia [--cache warm|cold] [--only r_seurat_bpcells,...]
#          [--include-optional] [--resume] [--run-id ID] [--single "feature|data|10:0|all" --rep N] [--start-time SECS]
#
# The SCUBA package is loaded from `[suite].r_scuba_path` (pkgload::load_all) so the development
# branch is benchmarked; SCUBA_JL_PATH is pointed at this checkout for the Julia engine.
T_PROCESS_START <- Sys.time()
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "harness.R"))
opts <- parse_cli(commandArgs(trailingOnly = TRUE))
ctx <- load_ctx(opts)
cfg <- ctx$cfg

## package setup ----------------------------------------------------------------------------------
t0 <- Sys.time()
scuba_path <- path.expand(cfg$suite$r_scuba_path)
suppressPackageStartupMessages(pkgload::load_all(scuba_path, quiet = TRUE))
Sys.setenv(SCUBA_JL_PATH = normalizePath(file.path(bench_dir(), "..")))
options(SCUBA.engine = opts$engine)
load_ms <- as.numeric(Sys.time() - t0, units = "secs") * 1000
suppressPackageStartupMessages({ library(SeuratObject); library(SingleCellExperiment); library(HDF5Array); library(BPCells) })
T_LOADED <- Sys.time()

STORAGE <- c(r_seurat_bpcells = "bpcells", r_seurat_bpcells_rowmajor = "bpcells_rowmajor", r_sce_backed = "sce_h5sparse",
             r_anndata_backed = "h5ad_csc", r_h5ad_julia = "h5ad_csc", r_seurat_inmem = "memory", raw_bpcells_r = "bpcells",
             raw_bpcells_r_rowmajor = "bpcells_rowmajor", raw_seurat_fetchdata = "memory")
# how each implementation names things
KEYS <- list(r_seurat_bpcells = list(prefix = "rna_", layers = c(data = "data", counts = "counts"), red = c(X_umap = "umap", X_pca = "pca")),
             r_seurat_bpcells_rowmajor = list(prefix = "rna_", layers = c(data = "data", counts = "counts"), red = c(X_umap = "umap", X_pca = "pca")),
             r_seurat_inmem = list(prefix = "rna_", layers = c(data = "data", counts = "counts"), red = c(X_umap = "umap", X_pca = "pca")),
             r_sce_backed = list(prefix = "RNA_", layers = c(data = "logcounts", counts = "counts"), red = c(X_umap = "UMAP", X_pca = "PCA")),
             r_anndata_backed = list(prefix = "X_", layers = c(data = NULL, counts = "counts"), red = c(X_umap = "X_umap", X_pca = "X_pca")),
             r_h5ad_julia = list(prefix = "X_", layers = c(data = "X", counts = "counts"), red = c(X_umap = "X_umap", X_pca = "X_pca")))

anndata_ok <- function() requireNamespace("anndata", quietly = TRUE) && packageVersion("anndata") < "0.8.0"

open_object <- function(impl) {
  d <- ctx$data
  switch(impl,
    r_seurat_bpcells = readRDS(file.path(d, "seurat_bpcells.rds")),
    r_seurat_bpcells_rowmajor = readRDS(file.path(d, "seurat_bpcells_rowmajor.rds")),
    raw_bpcells_r_rowmajor = list(data = open_matrix_dir(file.path(d, "bpcells", "RNA_data_rowmajor")),
                                  counts = open_matrix_dir(file.path(d, "bpcells", "RNA_counts_rowmajor"))),
    r_sce_backed = readRDS(file.path(d, "sce_h5sparse.rds")),
    r_anndata_backed = {
      reticulate::use_condaenv(cfg$suite$conda_env, required = FALSE)
      anndata::read_h5ad(file.path(d, "csc.h5ad"), backed = "r")
    },
    r_h5ad_julia = SCUBA::open_h5ad(file.path(d, "csc.h5ad")),    # file-backed scuba_source (Julia engine)
    raw_bpcells_r = list(data = open_matrix_dir(file.path(d, "bpcells", "RNA_data")),
                         counts = open_matrix_dir(file.path(d, "bpcells", "RNA_counts"))),
    r_seurat_inmem = , raw_seurat_fetchdata = {
      counts <- as(open_matrix_dir(file.path(d, "bpcells", "RNA_counts")), "dgCMatrix")
      dat <- as(open_matrix_dir(file.path(d, "bpcells", "RNA_data")), "dgCMatrix")
      obs <- as.data.frame(arrow::read_parquet(file.path(d, "obs.parquet"))); rownames(obs) <- obs$cell_id; obs$cell_id <- NULL
      obj <- CreateSeuratObject(counts = counts, meta.data = obs, assay = "RNA", min.cells = 0, min.features = 0)
      obj[["RNA"]]$data <- dat
      umap <- as.matrix(arrow::read_parquet(file.path(d, "obsm_X_umap.parquet"))); rownames(umap) <- colnames(obj); colnames(umap) <- c("UMAP_1", "UMAP_2")
      obj[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = "RNA")
      obj
    },
    stop("unknown implementation ", impl))
}

## one workload, once -------------------------------------------------------------------------------
run_work <- function(impl, obj, w, all_ids) {
  k <- KEYS[[impl]]
  if (w$workload == "feature") {
    genes <- ctx$selection$gene_sets[[w$gene_set]]
    cells <- cells_for_subset(ctx, all_ids, w$cell_subset)
    e <- ctx$expected$feature[[paste(w$layer, w$gene_set, w$cell_subset, sep = "|")]]
    if (impl %in% c("raw_bpcells_r", "raw_bpcells_r_rowmajor")) {
      m <- obj[[w$layer]]
      block <- if (is.null(cells)) as.matrix(m[genes, , drop = FALSE]) else as.matrix(m[genes, cells, drop = FALSE])
      return(list(chk = feature_checksum(block), n = ncol(block), exp = e$sum))
    }
    if (impl == "raw_seurat_fetchdata") {
      df <- SeuratObject::FetchData(obj, vars = paste0("rna_", genes), layer = k$layers[[w$layer]], cells = cells)
      return(list(chk = feature_checksum(df), n = nrow(df), exp = e$sum))
    }
    layer <- if (w$layer == "data" && impl == "r_anndata_backed") NULL else k$layers[[w$layer]]
    df <- fetch_data(obj, vars = paste0(k$prefix, genes), layer = layer, cells = cells)
    return(list(chk = feature_checksum(df), n = nrow(df), exp = e$sum))
  }
  if (w$workload == "feature_fast") {        # fetch_feature(): BPCells fast path on the legacy engine
    genes <- ctx$selection$gene_sets[[w$gene_set]]
    cells <- cells_for_subset(ctx, all_ids, w$cell_subset)
    e <- ctx$expected$feature[[paste(w$layer, w$gene_set, w$cell_subset, sep = "|")]]
    df <- fetch_feature(obj, features = genes, layer = k$layers[[w$layer]], cells = cells)
    return(list(chk = feature_checksum(df), n = nrow(df), exp = e$sum))
  }
  if (w$workload == "matrix_sparse") {
    genes <- ctx$selection$gene_sets[[w$gene_set]]
    e <- ctx$expected$feature[[paste("data", w$gene_set, "all", sep = "|")]]
    if (impl %in% c("raw_bpcells_r", "raw_bpcells_r_rowmajor")) { m <- as(obj$data[genes, , drop = FALSE], "dgCMatrix"); return(list(chk = feature_checksum(m), n = ncol(m), exp = e$sum)) }
    if (ctx$engine == "julia" && exists("fetch_matrix", asNamespace("SCUBA"))) {      # features x cells dgCMatrix straight from Julia
      m <- SCUBA::fetch_matrix(obj, features = genes, layer = k$layers[["data"]], sparse = TRUE)
      return(list(chk = feature_checksum(m), n = ncol(m), exp = e$sum))
    }
    df <- fetch_data(obj, vars = paste0(k$prefix, genes), layer = if (impl == "r_anndata_backed") NULL else k$layers[["data"]])
    return(list(chk = feature_checksum(df), n = nrow(df), exp = e$sum))
  }
  if (w$workload == "metadata") {
    if (impl %in% c("raw_bpcells_r", "raw_bpcells_r_rowmajor", "raw_seurat_fetchdata")) return(NULL)
    if (w$extra == "full") { df <- as.data.frame(fetch_metadata(obj, full_table = TRUE)); return(list(chk = metadata_checksum_full(df), n = nrow(df), exp = ctx$expected$metadata$full)) }
    v <- fetch_metadata(obj, vars = w$extra, return_class = "vector")
    return(list(chk = metadata_checksum(v), n = length(v), exp = ctx$expected$metadata[[w$extra]]))
  }
  if (w$workload == "reduction") {
    if (impl %in% c("raw_bpcells_r", "raw_bpcells_r_rowmajor", "raw_seurat_fetchdata")) return(NULL)
    p <- strsplit(w$extra, "|", fixed = TRUE)[[1]]; dims <- as.integer(strsplit(p[2], "-")[[1]])
    df <- fetch_reduction(obj, reduction = k$red[[p[1]]], dims = dims[1]:dims[2])   # all cells, as the Julia harness does
    return(list(chk = feature_checksum(as.data.frame(df)), n = nrow(df), exp = ctx$expected$reduction[[w$extra]]))
  }
  if (w$workload == "fetch_cells") {
    if (impl %in% c("raw_bpcells_r", "raw_bpcells_r_rowmajor", "raw_seurat_fetchdata")) return(NULL)
    p <- strsplit(w$extra, ":", fixed = TRUE)[[1]]
    ids <- fetch_cells(obj, meta_var = p[1], meta_levels = p[2])
    return(list(chk = fetch_cells_checksum(ids), n = length(ids), exp = ctx$expected$fetch_cells[[w$extra]]))
  }
  if (w$workload == "unique_values") {
    if (impl %in% c("raw_bpcells_r", "raw_bpcells_r_rowmajor", "raw_seurat_fetchdata")) return(NULL)
    v <- unique_values(obj, var = w$extra)
    return(list(chk = unique_values_checksum(v), n = length(v), exp = ctx$expected$unique_values[[w$extra]]))
  }
  if (w$workload == "julia_rtt") {
    for (i in 1:50) JuliaCall::julia_eval("1")
    return(list(chk = "", n = 50L, exp = NULL))
  }
  stop("unknown workload ", w$workload)
}

all_cells_of <- function(impl, obj) {
  if (impl %in% c("raw_bpcells_r", "raw_bpcells_r_rowmajor")) colnames(obj$data) else if (impl == "raw_seurat_fetchdata") colnames(obj) else get_all_cells(obj)
}

run_impl <- function(impl) {
  storage <- STORAGE[[impl]]
  reps <- cfg$reps; warm <- reps$warm; timeout_ms <- reps$timeout_s * 1000
  rules <- unlist(cfg$skip$rules)
  if (impl == "r_h5ad_julia" && (ctx$engine != "julia" || !exists("open_h5ad", asNamespace("SCUBA")))) {
    ctx <<- result_row(ctx, impl, storage, "open", phase = "open", status = "skipped", note = "file-backed scuba_source needs the Julia engine and SCUBA::open_h5ad()")
    return(invisible())
  }
  if (impl == "r_anndata_backed" && !anndata_ok()) {
    ctx <<- result_row(ctx, impl, storage, "open", phase = "open", status = "skipped",
                       note = sprintf("R anndata %s returns python-class objects; install remotes::install_version('anndata', '0.7.5.6')",
                                      tryCatch(as.character(packageVersion("anndata")), error = function(e) "missing")))
    return(invisible())
  }
  single <- nzchar(opts$single)
  tm <- timed(function() open_object(impl))
  obj <- tm$value
  if (!single) ctx <<- result_row(ctx, impl, storage, "open", phase = "open", wall_ms = tm$wall_ms, cpu_ms = tm$cpu_ms,
                                  open_ms = tm$wall_ms, alloc_mb = tm$alloc_mb, note = class(obj)[1])
  all_ids <- all_cells_of(impl, obj)
  works <- if (single) list(parse_single(opts$single)) else workloads(ctx)
  if (!single && impl %in% c("r_seurat_bpcells", "r_seurat_bpcells_rowmajor") && ctx$engine == "r") {      # legacy fetch_feature fast path
    extra <- Filter(function(w) w$workload == "feature" && w$cell_subset %in% c("all", "rand10"), works)
    works <- c(works, lapply(extra, function(w) { w$workload <- "feature_fast"; w }))
  }
  if (!single && ctx$engine == "julia" && impl == names(STORAGE)[1])
    works <- c(works, list(list(workload = "julia_rtt", layer = "", gene_set = "", n_features = 0L, cell_subset = "", extra = "50 x julia_eval('1')")))
  seen <- character()
  for (w in works) {
    label <- work_label(w)
    rule <- skip_rule_matches(rules, impl, ctx$scale, w$n_features, w$workload)
    if (!is.null(rule)) { ctx <<- result_row(ctx, impl, storage, label, w$n_features, gene_set_index(w), w$cell_subset, w$layer, status = "skipped", note = paste("rule:", rule)); next }
    family <- paste(w$workload, w$layer)
    nrep <- if (single) 1L else warm
    rep_start <- if (single) opts$rep else 0L
    budget_exceeded <- FALSE
    for (rep in rep_start:(rep_start + nrep - (if (single) 1L else 0L))) {
      phase <- if (rep == 0) { if (family %in% seen) "warmup" else "first_call" } else "fetch"
      if (!single && rep == 0) seen <- c(seen, family)
      key <- row_key(impl, label, w$n_features, gene_set_index(w), w$cell_subset, w$layer, ctx$cache_state, phase, rep)
      if (opts$resume && key %in% ctx$existing) next
      if (budget_exceeded) { ctx <<- result_row(ctx, impl, storage, label, w$n_features, gene_set_index(w), w$cell_subset, w$layer, phase = phase, rep = rep, status = "skipped_extrapolated", note = "first call exceeded timeout budget"); next }
      res <- tryCatch({
        tm <- timed(function() suppressWarnings(run_work(impl, obj, w, all_ids)))
        if (is.null(tm$value)) NULL else list(tm = tm, status = "ok", note = w$extra)
      }, error = function(e) list(tm = list(value = list(chk = "", n = 0L, exp = NULL), wall_ms = NA, cpu_ms = NA, alloc_mb = NA),
                                  status = if (grepl("two-length", conditionMessage(e))) "skipped" else "error",
                                  note = if (grepl("two-length", conditionMessage(e))) "legacy fetch_reduction supports exactly two dims" else substr(conditionMessage(e), 1, 300)))
      if (is.null(res)) break                                        # workload not applicable to this implementation
      if (res$status == "skipped") { ctx <<- result_row(ctx, impl, storage, label, w$n_features, gene_set_index(w), w$cell_subset, w$layer, phase = phase, rep = rep, status = "skipped", note = res$note); break }
      v <- res$tm$value
      ok <- if (res$status == "ok") checksum_ok(v$chk, v$exp) else ""
      ctx <<- result_row(ctx, impl, storage, label, w$n_features, gene_set_index(w), w$cell_subset, w$layer, n_cells_out = v$n,
                         phase = phase, rep = rep, wall_ms = res$tm$wall_ms, cpu_ms = res$tm$cpu_ms, alloc_mb = res$tm$alloc_mb,
                         result_checksum = if (is.numeric(v$chk)) format(v$chk, digits = 17) else as.character(v$chk),
                         checksum_ok = ok, status = res$status, note = res$note)
      if (res$status == "ok" && ok == "false") message("checksum mismatch: ", impl, " ", label, " ", w$gene_set, " ", w$cell_subset)
      if (res$status == "ok" && isTRUE(res$tm$wall_ms > timeout_ms)) budget_exceeded <- TRUE
    }
  }
  rm(obj); invisible(gc())
}

## main ------------------------------------------------------------------------------------------------
impls <- if (length(opts$only)) opts$only else unlist(cfg$implementations[[if (ctx$engine == "julia") "r_julia" else "r_legacy"]])
if (!length(opts$only) && !opts$include_optional) impls <- setdiff(impls, unlist(cfg$implementations$optional))
first <- TRUE
for (impl in impls) {
  ctx <- open_output(ctx, sprintf("r_%s_%s_%s_%s", impl, ctx$engine, ctx$scale, ctx$cache_state))
  if (first && !nzchar(opts$single)) {
    t_start <- if (is.na(opts$start_time)) as.numeric(T_PROCESS_START) else opts$start_time
    ctx <- result_row(ctx, impl, STORAGE[[impl]], "startup", phase = "startup", rep = 0L, wall_ms = (as.numeric(T_LOADED) - t_start) * 1000,
                      note = sprintf("process start -> packages loaded (load_all SCUBA %.0f ms)", load_ms))
    if (ctx$engine == "julia") {
      tm <- timed(function() scuba_julia_available(quiet = TRUE))
      ctx <- result_row(ctx, impl, STORAGE[[impl]], "startup", phase = "startup", rep = 1L, wall_ms = tm$wall_ms, cpu_ms = tm$cpu_ms,
                        note = paste("scuba_julia_available() JuliaCall init; engine =", scuba_engine()))
    }
    first <- FALSE
  }
  message(sprintf("[r %s] %s engine=%s scale=%s cache=%s", format(Sys.time(), "%H:%M:%S"), impl, ctx$engine, ctx$scale, ctx$cache_state))
  run_impl(impl)
  ctx <- close_output(ctx)
}
cat("results in", ctx$out_dir, "\n")
