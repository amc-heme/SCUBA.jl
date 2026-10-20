#!/usr/bin/env Rscript
# Generate golden outputs from the R package SCUBA (expected: v1.2.2) for every case in
# cases.json, on the Seurat and SCE forms of the AML fixtures. The AnnData form is skipped
# because the R anndata >= 0.8 package returns python objects that SCUBA 1.2.2 cannot dispatch on.
#
# Usage: Rscript fixtures/make_goldens.R [outdir]
suppressPackageStartupMessages({
  library(SCUBA); library(jsonlite); library(SingleCellExperiment)
})
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else file.path("fixtures", "build", "goldens")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cases <- fromJSON(file.path("fixtures", "cases.json"), simplifyVector = FALSE)
objects <- list(seurat = SCUBA::AML_Seurat, sce = SCUBA::AML_SCE())
cat("SCUBA version:", as.character(packageVersion("SCUBA")), "\n")

capture <- function(expr) {
  warnings <- character(); error <- NULL; value <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) { error <<- conditionMessage(e); NULL }),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(value = value, warnings = warnings, error = error)
}
col_types <- function(df) {
  lapply(df, function(x) {
    if (is.factor(x)) list(class = "factor", levels = levels(x), ordered = is.ordered(x))
    else list(class = class(x)[1])
  })
}
write_table <- function(df, path_base) {
  df <- as.data.frame(df)
  out <- cbind(cell = rownames(df), df, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  write.csv(out, paste0(path_base, ".csv"), row.names = FALSE, na = "")
  write_json(col_types(df), paste0(path_base, ".types.json"), auto_unbox = TRUE, pretty = TRUE)
}

conditions <- list()
for (cs in cases$cases) {
  for (form in cases$forms) {
    obj <- objects[[form]]
    cells <- NULL
    if (!is.null(cs$cells) && cs$cells == "first20_reversed") cells <- rev(get_all_cells(obj)[1:20])
    vars <- if (!is.null(cs$vars)) unlist(cs$vars[[form]]) else NULL
    layer <- if (!is.null(cs$layer)) cs$layer[[form]] else NULL
    res <- switch(cs$fn,
      fetch_data = capture(fetch_data(obj, vars = vars, layer = layer, cells = cells)),
      fetch_reduction = capture(fetch_reduction(obj, reduction = default_reduction(obj), cells = get_all_cells(obj))),
      fetch_metadata_full = capture(as.data.frame(fetch_metadata(obj, full_table = TRUE))),
      fetch_metadata = capture(fetch_metadata(obj, vars = vars)),
      fetch_cells = capture(fetch_cells(obj, meta_var = cs$meta_var, meta_levels = unlist(cs$meta_levels))),
      unique_values = capture(as.character(unique_values(obj, var = cs$var))),
      meta_varnames = capture(meta_varnames(obj)),
      stop("unknown fn ", cs$fn))
    base <- file.path(outdir, paste0(cs$id, "__", form))
    if (!is.null(res$value)) {
      if (is.data.frame(res$value)) write_table(res$value, base)
      else writeLines(as.character(res$value), paste0(base, ".txt"))
    }
    conditions[[paste0(cs$id, "__", form)]] <- list(warnings = as.list(res$warnings), error = res$error,
                                                     ncol = if (is.data.frame(res$value)) ncol(res$value) else NULL,
                                                     colnames = if (is.data.frame(res$value)) colnames(res$value) else NULL)
    cat(sprintf("%-32s %-7s %s\n", cs$id, form, if (is.null(res$error)) "ok" else paste("ERROR:", substr(res$error, 1, 60))))
  }
}
write_json(conditions, file.path(outdir, "conditions.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("wrote", length(conditions), "golden cases to", outdir, "\n")
