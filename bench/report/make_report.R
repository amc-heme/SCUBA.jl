#!/usr/bin/env Rscript
# Build figures and summary.md from results/results.parquet (run tools/merge_results.py first).
#   Rscript bench/report/make_report.R [--results bench/results/results.parquet] [--out bench/report]
suppressPackageStartupMessages({ library(arrow); library(ggplot2); library(dplyr); library(tidyr) })
args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
bench <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ".."))
res_path <- get("--results", file.path(bench, "results", "results.parquet"))
out <- get("--out", file.path(bench, "report")); fig <- file.path(out, "figures"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)

df <- read_parquet(res_path) %>% mutate(across(c(wall_ms, cpu_ms, open_ms, rss_mb, maxrss_mb, alloc_mb, n_features, rep), as.numeric))
ok <- df %>% filter(status == "ok", is.na(checksum_ok) | checksum_ok != FALSE | checksum_ok == "" | checksum_ok == TRUE)
ok <- df %>% filter(status == "ok", !(checksum_ok %in% c(FALSE)))
scale_cells <- c("10k" = 1e4, "100k" = 1e5, "1M" = 1e6, "3M" = 3e6)
ok$cells <- scale_cells[ok$scale]
ok$impl <- ifelse(ok$language == "r", paste0(ok$implementation, " [", ok$engine, "]"), ok$implementation)
ok$lang <- factor(c(julia = "Julia (SCUBA.jl)", r = "R (SCUBA)", python = "Python")[ok$language], levels = c("Julia (SCUBA.jl)", "R (SCUBA)", "Python"))
x_cells <- scale_x_log10(breaks = c(1e4, 1e5, 1e6, 3e6), labels = c("10k", "100k", "1M", "3M"))

# Fixed implementation -> colour slot assignment (reference palette, categorical slots in fixed order per language)
LANG_LABEL <- c(julia = "Julia (SCUBA.jl)", r = "R (SCUBA)", python = "Python")
PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
impl_levels <- c("jl_h5ad_csc", "jl_h5ad_csr", "jl_bpcells", "jl_bpcells_rowmajor", "jl_snapshot",
                 "r_seurat_bpcells [r]", "r_seurat_bpcells_rowmajor [r]", "r_sce_backed [r]", "raw_bpcells_r [r]", "raw_bpcells_r_rowmajor [r]",
                 "r_seurat_bpcells [julia]", "r_seurat_bpcells_rowmajor [julia]", "r_sce_backed [julia]",
                 "r_anndata_backed [r]", "r_seurat_inmem [r]", "raw_seurat_fetchdata [r]",
                 "raw_anndata_py", "raw_h5py_py", "py_scuba_sc")
slot <- c(1, 2, 3, 4, 5,  1, 2, 3, 4, 5,  6, 7, 8,  1, 2, 3,  1, 2, 3)
impl_cols <- setNames(PAL[slot], impl_levels)
present <- intersect(impl_levels, unique(ok$impl)); extra <- setdiff(unique(ok$impl), impl_levels)
ok$impl <- factor(ok$impl, levels = c(present, extra))
cols <- c(impl_cols[present], setNames(rep("#52514e", length(extra)), extra))

theme_bench <- function() theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.3),
        text = element_text(colour = "#0b0b0b"), axis.text = element_text(colour = "#52514e"), strip.text = element_text(face = "bold"),
        legend.position = "bottom", plot.title.position = "plot", plot.background = element_rect(fill = "#fcfcfb", colour = NA))
save <- function(p, name, w = 10, h = 6) ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 150, bg = "#fcfcfb")

summ <- ok %>% filter(phase == "fetch") %>%
  group_by(language, lang, impl, implementation, engine, scale, cells, workload, n_features, cell_subset, layer, cache_state) %>%
  summarise(median_ms = median(wall_ms, na.rm = TRUE), q25 = quantile(wall_ms, .25, na.rm = TRUE), q75 = quantile(wall_ms, .75, na.rm = TRUE),
            min_ms = min(wall_ms, na.rm = TRUE), n = n(), .groups = "drop")

## 1. Scaling curves: feature fetches (one figure per cell subset; panels = language x genes) --------
feat <- summ %>% filter(workload == "feature", layer == "data", cache_state == "warm", cell_subset %in% c("all", "rand10"))
if (nrow(feat) > 1 && length(unique(feat$scale)) > 1) {
  for (sub in c("all", "rand10")) {
    d <- feat %>% filter(cell_subset == sub)
    p <- ggplot(d, aes(cells, median_ms, colour = impl, group = impl)) +
      geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
      facet_grid(lang ~ n_features, labeller = labeller(n_features = function(x) paste(x, "genes"))) +
      x_cells + scale_y_log10(labels = scales::label_comma()) +
      scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
      labs(title = sprintf("Feature fetch time vs dataset size, %s (median of warm repetitions)", c(all = "all cells", rand10 = "10%% random cells")[sub]),
           x = "cells", y = "wall time (ms, log)") + theme_bench() + guides(colour = guide_legend(ncol = 4))
    save(p, paste0("scaling_features_", sub), 12, 8)
  }
}
## 2. Per-scale bars: feature fetch by n_features (all cells), panels per language -----------------------
for (sc in unique(feat$scale)) {
  d <- feat %>% filter(scale == sc, cell_subset == "all")
  p <- ggplot(d, aes(factor(n_features), median_ms, fill = impl)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    geom_errorbar(aes(ymin = q25, ymax = q75), position = position_dodge(width = 0.85), width = 0.25, linewidth = 0.3, colour = "#52514e") +
    facet_wrap(~lang, nrow = 1) +
    scale_y_log10(labels = scales::label_comma()) + scale_fill_manual(values = cols, name = NULL) +
    labs(title = sprintf("Feature fetch, all cells, %s cells (median, IQR)", sc), x = "genes requested", y = "wall time (ms, log)") + theme_bench() + guides(fill = guide_legend(ncol = 4))
  save(p, paste0("features_", sc), 12, 6)
}
## 3. Metadata / reductions / cells --------------------------------------------------------------------
other <- summ %>% filter(grepl("^(metadata|reduction|fetch_cells|unique_values)", workload), cache_state == "warm")
if (nrow(other)) {
  for (sc in unique(other$scale)) {
    d <- other %>% filter(scale == sc)
    p <- ggplot(d, aes(workload, median_ms, fill = impl)) + geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      facet_wrap(~lang, nrow = 1) +
      scale_y_log10(labels = scales::label_comma()) + scale_fill_manual(values = cols, name = NULL) + coord_flip() +
      labs(title = sprintf("Metadata, reductions and cell queries, %s cells (median)", sc), x = NULL, y = "wall time (ms, log)") + theme_bench() + guides(fill = guide_legend(ncol = 4))
    save(p, paste0("other_", sc), 11, 7)
  }
}
## 4. Cold vs warm ---------------------------------------------------------------------------------------
cw <- ok %>% filter(phase %in% c("fetch", "open"), cache_state %in% c("warm", "cold")) %>%
  group_by(impl, scale, workload, n_features, cell_subset, cache_state, phase) %>% summarise(median_ms = median(wall_ms, na.rm = TRUE), n = n(), .groups = "drop")
cold_pairs <- cw %>% filter(cache_state == "cold") %>% distinct(impl, scale, workload, n_features, cell_subset)
if (nrow(cold_pairs)) {
  d <- cw %>% semi_join(cold_pairs, by = c("impl", "scale", "workload", "n_features", "cell_subset")) %>% filter(phase == "fetch")
  p <- ggplot(d, aes(interaction(workload, n_features, cell_subset, sep = " "), median_ms, fill = cache_state)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) + facet_grid(scale ~ impl) + coord_flip() +
    scale_y_log10(labels = scales::label_comma()) + scale_fill_manual(values = c(warm = "#2a78d6", cold = "#eb6834"), name = NULL) +
    labs(title = "Cold (page cache evicted) vs warm fetch time (median)", x = NULL, y = "wall time (ms, log)") + theme_bench()
  save(p, "cold_vs_warm", 13, 7)
}
## 5. Memory ----------------------------------------------------------------------------------------------
mem <- ok %>% filter(phase == "fetch", cache_state == "warm", workload == "feature", n_features == 1000, cell_subset == "all") %>%
  group_by(lang, impl, scale) %>% summarise(maxrss_mb = max(maxrss_mb, na.rm = TRUE), alloc_mb = median(alloc_mb, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(c(maxrss_mb, alloc_mb), names_to = "metric", values_to = "mb") %>% filter(is.finite(mb))
if (nrow(mem)) {
  p <- ggplot(mem, aes(scale, mb, fill = impl)) + geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    facet_grid(metric ~ lang, scales = "free_y", labeller = labeller(metric = c(maxrss_mb = "peak process RSS (MB)", alloc_mb = "allocation during a 1000-gene fetch (MB)"))) +
    scale_fill_manual(values = cols, name = NULL) + labs(title = "Memory", x = "cells", y = "MB") + theme_bench() + guides(fill = guide_legend(ncol = 4))
  save(p, "memory", 12, 6)
}
## 6. Concurrency --------------------------------------------------------------------------------------------
conc_files <- list.files(file.path(bench, "results", "raw"), pattern = "concurrency_summary.csv", recursive = TRUE, full.names = TRUE)
if (length(conc_files)) {
  cs <- bind_rows(lapply(conc_files, function(f) { d <- read.csv(f); if (!"scale" %in% names(d)) d$scale <- "unknown"; if (!"engine" %in% names(d)) d$engine <- "native"; d })) %>%
    filter(workload == "all") %>%
    mutate(impl = paste0(implementation, ifelse(language == "r", paste0(" [", engine, "]"), ""), ifelse(mode == "threads", paste0(" threads/", share), "")),
           lang = factor(c(julia = "Julia (SCUBA.jl)", r = "R (SCUBA)", python = "Python")[language], levels = c("Julia (SCUBA.jl)", "R (SCUBA)", "Python")))
  conc_cols <- setNames(PAL[((seq_along(unique(cs$impl)) - 1) %% 8) + 1], unique(cs$impl))
  p1 <- ggplot(cs, aes(n, throughput_qps, colour = impl, group = impl)) + geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    facet_grid(scale ~ lang, scales = "free_y") + scale_y_log10(labels = scales::label_comma()) +
    scale_x_continuous(breaks = unique(cs$n)) + labs(title = "Throughput vs concurrent sessions (mixed query stream)", x = "concurrent sessions", y = "queries / s (log)") + theme_bench() + scale_colour_manual(values = conc_cols, name = NULL) + guides(colour = guide_legend(ncol = 3))
  p2 <- ggplot(cs, aes(n, p95_ms, colour = impl, group = impl)) + geom_line(linewidth = 0.7) + geom_point(size = 2.2) + scale_y_log10(labels = scales::label_comma()) +
    facet_grid(scale ~ lang) +
    scale_x_continuous(breaks = unique(cs$n)) + labs(title = "p95 latency vs concurrent sessions", x = "concurrent sessions", y = "p95 wall time (ms, log)") + theme_bench() + scale_colour_manual(values = conc_cols, name = NULL) + guides(colour = guide_legend(ncol = 3))
  save(p1, "concurrency_throughput", 11, 7); save(p2, "concurrency_p95", 11, 7)
}

## 7. summary.md -----------------------------------------------------------------------------------------------
md <- c("# SCUBA benchmark summary", "", sprintf("Generated %s from `%s` (%d rows, %d ok).", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), basename(res_path), nrow(df), nrow(ok)), "")
env_files <- list.files(file.path(bench, "results", "env"), full.names = TRUE)
if (length(env_files)) {
  e <- jsonlite::fromJSON(tail(sort(env_files), 1))
  md <- c(md, "## Environment", "", sprintf("- host: %s (%s)", e$host, e$platform),
          sprintf("- hardware: %s, %s cores, %.0f GB RAM", e$hardware$hw.model, e$hardware$hw.ncpu, as.numeric(e$hardware$hw.memsize) / 1e9),
          sprintf("- julia: %s; SCUBA.jl %s%s", e$julia$version, substr(e$scuba_jl$sha, 1, 8), if (isTRUE(e$scuba_jl$dirty)) " (dirty)" else ""),
          sprintf("- R: %s; SCUBA R %s (%s)%s", e$R$version, substr(e$scuba_r$sha, 1, 8), e$scuba_r$branch, if (isTRUE(e$scuba_r$dirty)) " (dirty)" else ""),
          sprintf("- python env `%s`: %s", e$python$env, paste(names(e$python$packages), unlist(e$python$packages), sep = " ", collapse = ", ")), "")
}
fmt_ms <- function(x) ifelse(is.na(x), "", ifelse(x >= 1000, sprintf("%.1f s", x / 1000), sprintf("%.1f ms", x)))
for (sc in intersect(names(scale_cells), unique(summ$scale))) {
  d <- summ %>% filter(scale == sc, cache_state == "warm", workload == "feature", layer == "data", cell_subset == "all") %>%
    mutate(cell = sprintf("%s (IQR %s–%s)", fmt_ms(median_ms), fmt_ms(q25), fmt_ms(q75))) %>% select(impl, n_features, cell) %>%
    pivot_wider(names_from = n_features, values_from = cell, names_prefix = "genes=")
  md <- c(md, sprintf("## %s cells: feature fetch over all cells, warm (median, IQR)", sc), "", paste0("| implementation | ", paste(setdiff(names(d), "impl"), collapse = " | "), " |"),
          paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
          apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
  o <- summ %>% filter(scale == sc, cache_state == "warm", grepl("^(metadata|reduction|fetch_cells|unique_values)", workload)) %>%
    mutate(cell = fmt_ms(median_ms)) %>% select(impl, workload, cell) %>% pivot_wider(names_from = workload, values_from = cell)
  if (nrow(o)) md <- c(md, sprintf("## %s cells: metadata, reductions, cells (median)", sc), "", paste0("| implementation | ", paste(setdiff(names(o), "impl"), collapse = " | "), " |"),
                       paste0("|", paste(rep("---", ncol(o)), collapse = "|"), "|"), apply(o, 1, function(r) paste0("| ", paste(ifelse(is.na(r), "", r), collapse = " | "), " |")), "")
}
st <- ok %>% filter(phase %in% c("startup", "open", "first_call")) %>% group_by(impl, scale, phase) %>% summarise(ms = median(wall_ms, na.rm = TRUE), .groups = "drop") %>%
  mutate(ms = fmt_ms(ms)) %>% pivot_wider(names_from = phase, values_from = ms)
md <- c(md, "## Startup, open and first-call latency (median)", "", paste0("| implementation | scale | ", paste(setdiff(names(st), c("impl", "scale")), collapse = " | "), " |"),
        paste0("|", paste(rep("---", ncol(st)), collapse = "|"), "|"), apply(st, 1, function(r) paste0("| ", paste(ifelse(is.na(r), "", r), collapse = " | "), " |")), "")
sk <- df %>% filter(status != "ok") %>% count(implementation, engine, scale, status, note) %>% arrange(desc(n))
if (nrow(sk)) md <- c(md, "## Skipped and failed cells", "", "| implementation | engine | scale | status | note | rows |", "|---|---|---|---|---|---|",
                      apply(head(sk, 60), 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
bad <- df %>% filter(status == "ok", checksum_ok %in% c(FALSE)) %>% count(implementation, engine, scale, workload)
md <- c(md, "## Checksum mismatches", "", if (nrow(bad)) c("| implementation | engine | scale | workload | rows |", "|---|---|---|---|---|", apply(bad, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))) else "none", "")
md <- c(md, "## Figures", "", paste0("- `figures/", list.files(fig, pattern = "png$"), "`"))
writeLines(md, file.path(out, "summary.md"))
cat("wrote", file.path(out, "summary.md"), "and", length(list.files(fig, pattern = "png$")), "figures\n")
