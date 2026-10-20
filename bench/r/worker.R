#!/usr/bin/env Rscript
# Concurrency worker (R): one simulated Shiny session issuing a mixed query stream.
#   Rscript bench/r/worker.R --scale 1M --impl r_seurat_bpcells --engine r|julia --worker 2 --duration 60 --go <barrier file> [--run-id ID] [--sessions N]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "harness.R"))
args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
opts <- list(scale = get("--scale"), config = get("--config", file.path(bench_dir(), "config", "bench.toml")), cache = "warm",
             only = character(), include_optional = FALSE, resume = FALSE, run_id = get("--run-id", ""), single = "", rep = 1L,
             engine = get("--engine", "r"), start_time = NA_real_)
impl <- get("--impl"); worker <- as.integer(get("--worker", "1")); duration <- as.numeric(get("--duration", "60"))
go <- get("--go", ""); seed <- as.integer(get("--seed", "7")); n_workers <- as.integer(get("--sessions", "1")); warmup <- as.numeric(get("--warmup", "10"))
ctx <- load_ctx(opts); cfg <- ctx$cfg
suppressPackageStartupMessages(pkgload::load_all(path.expand(cfg$suite$r_scuba_path), quiet = TRUE))
Sys.setenv(SCUBA_JL_PATH = normalizePath(file.path(bench_dir(), "..")))
options(SCUBA.engine = opts$engine)
suppressPackageStartupMessages({ library(SeuratObject); library(SingleCellExperiment); library(HDF5Array); library(BPCells) })
if (opts$engine == "julia") invisible(scuba_julia_available(quiet = TRUE))
# reuse the harness run.R definitions without executing its main block
src_lines <- readLines(file.path(bench_dir(), "r", "run.R"))
main_at <- grep("^## main", src_lines)
eval(parse(text = src_lines[(grep("^STORAGE <- ", src_lines)[1]):(main_at - 1)]))
ctx <- open_output(ctx, sprintf("r_conc_%s_%s_%s_process_n%d_w%d", impl, opts$engine, ctx$scale, n_workers, worker))
obj <- open_object(impl); all_ids <- all_cells_of(impl, obj)
set.seed(seed + worker)
mix <- cfg$concurrency$mix; mnames <- names(mix); mw <- unlist(mix)
pick <- function() {
  name <- sample(mnames, 1, prob = mw)
  switch(name,
    feature1_all = list(workload = "feature", layer = "data", gene_set = sprintf("1:%d", sample(0:2, 1)), n_features = 1L, cell_subset = "all", extra = ""),
    feature10_rand10 = list(workload = "feature", layer = "data", gene_set = sprintf("10:%d", sample(0:2, 1)), n_features = 10L, cell_subset = "rand10", extra = ""),
    meta_col = list(workload = "metadata", layer = "", gene_set = "", n_features = 0L, cell_subset = "all", extra = sample(c("cell_type", "nCount_RNA"), 1)),
    reduction = list(workload = "reduction", layer = "", gene_set = "", n_features = 0L, cell_subset = "all", extra = "X_umap|1-2"),
    fetch_cells = list(workload = "fetch_cells", layer = "", gene_set = "", n_features = 0L, cell_subset = "", extra = sample(unlist(cfg$workloads$fetch_cells), 1)))
}
t_warm <- Sys.time()
while (if (nzchar(go)) !file.exists(go) else as.numeric(Sys.time() - t_warm, units = "secs") < warmup) {
  invisible(suppressWarnings(run_work(impl, obj, pick(), all_ids)))
  if (as.numeric(Sys.time() - t_warm, units = "secs") > 300) break
}
t_end <- Sys.time() + duration; n <- 0L
while (Sys.time() < t_end) {
  w <- pick(); tm <- timed(function() suppressWarnings(run_work(impl, obj, w, all_ids))); v <- tm$value
  ok <- if (runif(1) < cfg$concurrency$checksum_fraction) checksum_ok(v$chk, v$exp) else ""
  ctx <- result_row(ctx, impl, STORAGE[[impl]], work_label(w), w$n_features, gene_set_index(w), w$cell_subset, w$layer, n_cells_out = v$n,
                    phase = "concurrent", rep = worker, wall_ms = tm$wall_ms, cpu_ms = tm$cpu_ms, alloc_mb = tm$alloc_mb,
                    result_checksum = if (is.numeric(v$chk)) format(v$chk, digits = 17) else as.character(v$chk), checksum_ok = ok,
                    note = sprintf("mode=process;n=%d", n_workers), threads = n_workers)
  n <- n + 1L
}
ctx <- close_output(ctx)
cat(sprintf("worker %d: %d queries\n", worker, n))
