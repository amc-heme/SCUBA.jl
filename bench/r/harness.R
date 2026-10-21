# Shared helpers for the R benchmark harness (schema in common/schema.md).
suppressPackageStartupMessages({ library(RcppTOML); library(jsonlite) })

SCHEMA <- c("suite_version", "run_id", "timestamp", "host", "language", "implementation", "engine", "storage",
            "scale", "workload", "n_features", "gene_set", "cell_subset", "layer", "n_cells_out", "cache_state",
            "phase", "rep", "threads", "wall_ms", "cpu_ms", "open_ms", "rss_mb", "maxrss_mb", "alloc_mb", "gc_ms",
            "result_checksum", "checksum_ok", "status", "note")

bench_dir <- function() {
  f <- grep("--file=", commandArgs(), value = TRUE)
  if (length(f)) normalizePath(file.path(dirname(sub("--file=", "", f[1])), "..")) else normalizePath(".")
}

parse_cli <- function(args) {
  opts <- list(config = file.path(bench_dir(), "config", "bench.toml"), cache = "warm", only = character(),
               include_optional = FALSE, resume = FALSE, run_id = "", single = "", rep = 1L, engine = "r",
               start_time = NA_real_)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--include-optional") opts$include_optional <- TRUE
    else if (a == "--resume") opts$resume <- TRUE
    else {
      v <- args[i + 1]; i <- i + 1
      switch(a,
        "--scale" = opts$scale <- v, "--config" = opts$config <- v, "--cache" = opts$cache <- v,
        "--only" = opts$only <- strsplit(v, ",")[[1]], "--run-id" = opts$run_id <- v, "--single" = opts$single <- v,
        "--rep" = opts$rep <- as.integer(v), "--engine" = opts$engine <- v, "--start-time" = opts$start_time <- as.numeric(v),
        stop("unknown argument ", a))
    }
    i <- i + 1
  }
  if (is.null(opts$scale)) stop("--scale is required")
  opts
}

load_ctx <- function(opts) {
  cfg <- parseTOML(opts$config)
  d <- file.path(path.expand(cfg$suite$data_dir), opts$scale)
  man <- fromJSON(file.path(d, "MANIFEST.json"))
  stopifnot(isTRUE(man$complete))
  sel <- fromJSON(file.path(d, "selection.json"), simplifyVector = TRUE)
  exp <- fromJSON(file.path(d, "expected.json"), simplifyVector = FALSE)
  run_id <- if (nzchar(opts$run_id)) opts$run_id else if (nzchar(cfg$suite$run_id)) cfg$suite$run_id else format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
  out_dir <- file.path(bench_dir(), cfg$suite$results_dir, "raw", run_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  list(cfg = cfg, scale = opts$scale, run_id = run_id, data = d, selection = sel, expected = exp, manifest = man,
       cache_state = opts$cache, out_dir = out_dir, host = Sys.info()[["nodename"]], con = NULL, existing = character(),
       engine = opts$engine, peak_rss = 0)
}

row_key <- function(...) paste(..., sep = "|")

existing_keys <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(character())
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE, colClasses = "character"), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(character())
  with(df, row_key(implementation, workload, n_features, gene_set, cell_subset, layer, cache_state, phase, rep))
}

open_output <- function(ctx, name) {
  path <- file.path(ctx$out_dir, paste0(name, ".csv"))
  fresh <- !file.exists(path) || file.size(path) == 0
  ctx$existing <- if (fresh) character() else existing_keys(path)
  ctx$con <- file(path, open = "a")
  if (fresh) writeLines(paste(SCHEMA, collapse = ","), ctx$con)
  flush(ctx$con)
  ctx$out_path <- path
  ctx
}
close_output <- function(ctx) { if (!is.null(ctx$con)) close(ctx$con); ctx$con <- NULL; ctx }

csv_field <- function(x) {
  s <- if (is.null(x) || length(x) == 0 || is.na(x)) "" else as.character(x)
  if (grepl('[,"\n]', s)) paste0('"', gsub('"', '""', s), '"') else s
}
fmt <- function(x) if (is.numeric(x)) { if (is.na(x)) "" else format(round(x, 4), scientific = FALSE, trim = TRUE) } else x

rss_mb <- function() tryCatch(as.numeric(ps::ps_memory_info()[["rss"]]) / 2^20, error = function(e) NA_real_)

result_row <- function(ctx, implementation, storage, workload, n_features = 0L, gene_set = -1L, cell_subset = "",
                       layer = "", n_cells_out = 0L, phase = "fetch", rep = 0L, wall_ms = NA_real_, cpu_ms = NA_real_,
                       open_ms = 0, alloc_mb = NA_real_, gc_ms = NA_real_, result_checksum = "", checksum_ok = "",
                       status = "ok", note = "", threads = 1L) {
  rss <- rss_mb(); ctx$peak_rss <- max(ctx$peak_rss, rss, na.rm = TRUE)
  vals <- list(suite_version = ctx$cfg$suite$version, run_id = ctx$run_id,
               timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), host = ctx$host, language = "r",
               implementation = implementation, engine = ctx$engine, storage = storage, scale = ctx$scale,
               workload = workload, n_features = n_features, gene_set = gene_set, cell_subset = cell_subset,
               layer = layer, n_cells_out = n_cells_out, cache_state = ctx$cache_state, phase = phase, rep = rep,
               threads = threads, wall_ms = fmt(wall_ms), cpu_ms = fmt(cpu_ms), open_ms = fmt(open_ms), rss_mb = fmt(rss),
               maxrss_mb = fmt(ctx$peak_rss), alloc_mb = fmt(alloc_mb), gc_ms = fmt(gc_ms),
               result_checksum = result_checksum, checksum_ok = checksum_ok, status = status, note = note)
  writeLines(paste(vapply(SCHEMA, function(c) csv_field(vals[[c]]), ""), collapse = ","), ctx$con)
  flush(ctx$con)
  ctx$existing <- c(ctx$existing, row_key(implementation, workload, n_features, gene_set, cell_subset, layer, ctx$cache_state, phase, rep))
  ctx
}

# Run f() once: list(value, wall_ms, cpu_ms, alloc_mb)
timed <- function(f) {
  invisible(gc(reset = TRUE, full = FALSE))
  t0 <- proc.time()
  value <- f()
  t1 <- proc.time()
  g <- gc(full = FALSE)
  alloc <- sum(g[, 6])   # "max used (Mb)" since the reset, Ncells + Vcells
  list(value = value, wall_ms = unname((t1 - t0)[["elapsed"]] * 1000),
       cpu_ms = unname(((t1 - t0)[["user.self"]] + (t1 - t0)[["sys.self"]]) * 1000), alloc_mb = alloc)
}

sha1_lines <- function(v) digest::digest(paste(as.character(v), collapse = "\n"), algo = "sha1", serialize = FALSE)

feature_checksum <- function(x) {
  if (is.data.frame(x)) sum(vapply(x, function(col) sum(as.numeric(col)), numeric(1)))
  else if (inherits(x, "dgCMatrix")) sum(as.numeric(x@x))
  else sum(as.numeric(x))
}
metadata_checksum_full <- function(df) sprintf("%d:%s", nrow(df), format(sum(as.numeric(df$nCount_RNA)), digits = 17))
metadata_checksum <- function(col) {
  if (is.numeric(col)) sprintf("%d:%s", length(col), format(sum(as.numeric(col), na.rm = TRUE), digits = 17))
  else sprintf("%d:%s", length(col), sha1_lines(as.character(col)))
}
fetch_cells_checksum <- function(ids) sprintf("%d:%s", length(ids), sha1_lines(sort(as.character(ids), method = "radix")))   # C-locale order like Python/Julia
unique_values_checksum <- function(v) sprintf("%d:%s", length(v), sha1_lines(sort(as.character(v), method = "radix")))

checksum_ok <- function(got, expected, rtol = 1e-5) {
  if (is.null(expected)) return("")
  if (is.numeric(got) && is.numeric(expected)) return(tolower(isTRUE(abs(got - expected) <= rtol * abs(expected) + 1e-9)))
  a <- strsplit(as.character(got), ":", fixed = TRUE)[[1]]; b <- strsplit(as.character(expected), ":", fixed = TRUE)[[1]]
  if (length(a) >= 2 && length(b) >= 2) {
    if (a[1] != b[1]) return("false")
    ga <- suppressWarnings(as.numeric(a[2])); gb <- suppressWarnings(as.numeric(b[2]))
    if (!is.na(ga) && !is.na(gb)) return(tolower(isTRUE(abs(ga - gb) <= rtol * abs(gb) + 1e-9)))
    return(tolower(a[2] == b[2]))
  }
  tolower(identical(as.character(got), as.character(expected)))
}

skip_rule_matches <- function(rules, impl, scale, n_features, workload) {
  for (r in rules) {
    m <- TRUE
    for (p in strsplit(r, " ")[[1]]) {
      if (startsWith(p, "impl=")) m <- m && sub("impl=", "", p) == impl
      else if (startsWith(p, "scale=")) m <- m && sub("scale=", "", p) == scale
      else if (startsWith(p, "n_features>=")) m <- m && n_features >= as.integer(sub("n_features>=", "", p))
      else if (startsWith(p, "workload=")) m <- m && sub("workload=", "", p) == workload
    }
    if (m) return(r)
  }
  NULL
}

cells_for_subset <- function(ctx, all_ids, sub) {
  if (sub == "all") return(NULL)
  s <- ctx$selection$cell_subsets[[sub]]
  if (startsWith(sub, "contig")) all_ids[(s[1] + 1):s[2]] else all_ids[s + 1]
}

# Workload list identical to the Julia harness
workloads <- function(ctx) {
  wl <- ctx$cfg$workloads
  W <- list()
  keys <- names(ctx$selection$gene_sets)
  ord <- order(as.integer(sub(":.*", "", keys)), as.integer(sub(".*:", "", keys)))
  for (key in keys[ord]) {
    n <- as.integer(sub(":.*", "", key)); rep <- as.integer(sub(".*:", "", key))
    layers <- unlist(wl$layers_all_sets); if (rep == 0) layers <- c(layers, unlist(wl$layers_first_set))
    for (layer in layers) for (sub in unlist(wl$cell_subsets))
      W[[length(W) + 1]] <- list(workload = "feature", layer = layer, gene_set = key, n_features = n, cell_subset = sub, extra = "")
    if (n == max(unlist(wl$n_features)))
      W[[length(W) + 1]] <- list(workload = "matrix_sparse", layer = "data", gene_set = key, n_features = n, cell_subset = "all", extra = "")
  }
  for (m in unlist(wl$metadata)) W[[length(W) + 1]] <- list(workload = "metadata", layer = "", gene_set = "", n_features = 0L, cell_subset = "all", extra = m)
  for (r in wl$reductions) W[[length(W) + 1]] <- list(workload = "reduction", layer = "", gene_set = "", n_features = 0L, cell_subset = "all",
                                                       extra = sprintf("%s|%d-%d", r$name, r$dims[[1]], r$dims[[2]]))
  for (f in unlist(wl$fetch_cells)) W[[length(W) + 1]] <- list(workload = "fetch_cells", layer = "", gene_set = "", n_features = 0L, cell_subset = "", extra = f)
  for (u in unlist(wl$unique_values)) W[[length(W) + 1]] <- list(workload = "unique_values", layer = "", gene_set = "", n_features = 0L, cell_subset = "", extra = u)
  W
}

parse_single <- function(spec) {
  p <- strsplit(spec, "|", fixed = TRUE)[[1]]
  kind <- p[1]
  if (kind %in% c("feature", "matrix_sparse")) list(workload = kind, layer = p[2], gene_set = p[3], n_features = as.integer(sub(":.*", "", p[3])), cell_subset = p[4], extra = "")
  else if (kind %in% c("metadata", "reduction")) list(workload = kind, layer = "", gene_set = "", n_features = 0L, cell_subset = "all", extra = paste(p[-1], collapse = "|"))
  else list(workload = kind, layer = "", gene_set = "", n_features = 0L, cell_subset = "", extra = p[2])
}
work_label <- function(w) if (w$workload %in% c("feature", "matrix_sparse")) w$workload else paste0(w$workload, ":", w$extra)
gene_set_index <- function(w) if (nzchar(w$gene_set)) as.integer(sub(".*:", "", w$gene_set)) else -1L
