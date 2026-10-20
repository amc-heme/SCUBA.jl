#!/usr/bin/env Rscript
# Build the R-side forms of a benchmark dataset from the generated h5ad files:
#   <data_dir>/<scale>/bpcells/RNA_counts, RNA_data   BPCells matrix directories (plain MatrixDir, packed)
#   <data_dir>/<scale>/seurat_bpcells.rds               Seurat v5 object whose RNA layers are the directories above
#   <data_dir>/<scale>/sce_h5sparse.rds                 SingleCellExperiment over HDF5Array::H5SparseMatrix(csc.h5ad)
# Usage: Rscript generate/make_r_objects.R --scale 10k [--config config/bench.toml] [--data-dir DIR] [--force]
suppressPackageStartupMessages({
  library(RcppTOML); library(arrow); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
scale <- get_arg("--scale"); stopifnot(!is.null(scale))
bench_dir <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ".."))
config <- get_arg("--config", file.path(bench_dir, "config", "bench.toml"))
cfg <- parseTOML(config)
data_root <- path.expand(get_arg("--data-dir", cfg$suite$data_dir))
force <- "--force" %in% args
d <- file.path(data_root, scale)
manifest <- fromJSON(file.path(d, "MANIFEST.json"))
stopifnot(isTRUE(manifest$complete))
log <- function(...) cat(sprintf("[make_r_objects %s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
timings <- list()

obs <- as.data.frame(read_parquet(file.path(d, "obs.parquet")))
cells <- obs$cell_id; obs$cell_id <- NULL; rownames(obs) <- cells
for (col in names(obs)) if (is.character(obs[[col]]) && col != "barcode") obs[[col]] <- factor(obs[[col]])
umap <- as.matrix(read_parquet(file.path(d, "obsm_X_umap.parquet"))); rownames(umap) <- cells
pca <- as.matrix(read_parquet(file.path(d, "obsm_X_pca.parquet"))); rownames(pca) <- cells
genes <- NULL

## 1. BPCells directories ---------------------------------------------------------------------
suppressPackageStartupMessages(library(BPCells))
src <- file.path(d, paste0(cfg$generation$bpcells_from, ".h5ad"))
bp_dir <- file.path(d, "bpcells")
counts_dir <- file.path(bp_dir, "RNA_counts"); data_dir <- file.path(bp_dir, "RNA_data")
if (force || !file.exists(file.path(counts_dir, "version"))) {
  t0 <- Sys.time()
  counts <- open_matrix_anndata_hdf5(src, group = "layers/counts")
  counts <- convert_matrix_type(counts, "uint32_t")
  write_matrix_dir(counts, counts_dir, compress = TRUE, overwrite = TRUE)
  dat <- open_matrix_anndata_hdf5(src, group = "X")
  write_matrix_dir(dat, data_dir, compress = TRUE, overwrite = TRUE)
  timings$bpcells_s <- as.numeric(Sys.time() - t0, units = "secs")
  log("wrote BPCells dirs in ", round(timings$bpcells_s, 1), "s; storage_order = ",
      readLines(file.path(counts_dir, "storage_order")), ", version = ", readLines(file.path(counts_dir, "version")))
} else log("BPCells dirs exist, skipping")
mcounts <- open_matrix_dir(counts_dir); mdata <- open_matrix_dir(data_dir)
genes <- rownames(mcounts)
stopifnot(identical(colnames(mcounts), cells), length(genes) == manifest$genes)
if (isTRUE(cfg$generation$bpcells_gene_major)) {
  gm <- file.path(bp_dir, "RNA_data_genemajor")
  if (force || !file.exists(file.path(gm, "version"))) {
    write_matrix_dir(transpose_storage_order(mdata), gm, compress = TRUE, overwrite = TRUE)
    log("wrote gene-major BPCells dir")
  }
}

## 2. Seurat v5 object with BPCells layers -----------------------------------------------------
seurat_path <- file.path(d, "seurat_bpcells.rds")
if (force || !file.exists(seurat_path)) {
  suppressPackageStartupMessages(library(SeuratObject))
  t0 <- Sys.time()
  obj <- CreateSeuratObject(counts = mcounts, meta.data = obs, assay = "RNA", min.cells = 0, min.features = 0)
  obj[["RNA"]]$data <- mdata
  u <- umap; colnames(u) <- c("UMAP_1", "UMAP_2")
  obj[["umap"]] <- CreateDimReducObject(embeddings = u, key = "UMAP_", assay = "RNA")
  p <- pca; colnames(p) <- paste0("PC_", seq_len(ncol(p)))
  obj[["pca"]] <- CreateDimReducObject(embeddings = p, key = "PC_", assay = "RNA")
  saveRDS(obj, seurat_path)
  timings$seurat_s <- as.numeric(Sys.time() - t0, units = "secs")
  log("wrote seurat_bpcells.rds in ", round(timings$seurat_s, 1), "s")
  rm(obj)
} else log("seurat_bpcells.rds exists, skipping")

## 3. SingleCellExperiment over HDF5-backed sparse assays --------------------------------------
sce_path <- file.path(d, "sce_h5sparse.rds")
if (force || !file.exists(sce_path)) {
  suppressPackageStartupMessages({ library(HDF5Array); library(SingleCellExperiment) })
  t0 <- Sys.time()
  csc <- file.path(d, "csc.h5ad")
  X <- H5SparseMatrix(csc, "X")
  C <- H5SparseMatrix(csc, "layers/counts")
  # HDF5Array presents an anndata csc_matrix group as genes x cells; older versions give cells x genes
  orient <- function(M) {
    if (identical(dim(M), c(length(genes), length(cells)))) M
    else if (identical(dim(M), c(length(cells), length(genes)))) t(M)
    else stop("unexpected H5SparseMatrix dims: ", paste(dim(M), collapse = "x"))
  }
  logcounts <- orient(X); counts_assay <- orient(C)
  dimnames(logcounts) <- list(genes, cells); dimnames(counts_assay) <- list(genes, cells)
  sce <- SingleCellExperiment(assays = list(counts = counts_assay, logcounts = logcounts),
                              colData = S4Vectors::DataFrame(obs),
                              reducedDims = list(PCA = pca, UMAP = umap), mainExpName = "RNA")
  base::saveRDS(sce, sce_path)   # BiocGenerics::saveRDS refuses out-of-memory seeds; the h5ad path is absolute
  timings$sce_s <- as.numeric(Sys.time() - t0, units = "secs")
  log("wrote sce_h5sparse.rds in ", round(timings$sce_s, 1), "s")
} else log("sce_h5sparse.rds exists, skipping")

## 4. Record in the manifest ---------------------------------------------------------------------
manifest$r_objects <- list(
  bpcells = list(counts = normalizePath(counts_dir), data = normalizePath(data_dir),
                 storage_order = readLines(file.path(counts_dir, "storage_order")),
                 version_counts = readLines(file.path(counts_dir, "version")), version_data = readLines(file.path(data_dir, "version"))),
  seurat = seurat_path, sce = sce_path, timings_s = timings,
  R = R.version.string, BPCells = as.character(packageVersion("BPCells")),
  SeuratObject = as.character(packageVersion("SeuratObject")), HDF5Array = as.character(packageVersion("HDF5Array")),
  created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
write_json(manifest, file.path(d, "MANIFEST.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
log("done")
