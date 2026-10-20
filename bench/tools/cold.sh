#!/usr/bin/env bash
# Evict a dataset's files from the page cache and verify.
#   tools/cold.sh <data_dir>/<scale> <storage>      storage: h5ad_csc | h5ad_csr | bpcells | seurat | sce | all
# Exit 0 when < 1% of pages remain resident, 2 when eviction could not be verified, 1 on usage error.
set -u
d="$1"; storage="${2:-all}"
case "$storage" in
  h5ad_csc) files=("$d/csc.h5ad") ;;
  h5ad_csr) files=("$d/csr.h5ad") ;;
  bpcells)  files=("$d/bpcells") ;;
  seurat)   files=("$d/seurat_bpcells.rds" "$d/bpcells") ;;
  sce)      files=("$d/sce_h5sparse.rds" "$d/csc.h5ad") ;;
  all)      files=("$d") ;;
  *) echo "unknown storage $storage" >&2; exit 1 ;;
esac
if sudo -n true 2>/dev/null; then sudo -n purge && echo "purged (sudo)"; fi
if ! command -v vmtouch >/dev/null; then echo "vmtouch not found (brew install vmtouch)" >&2; exit 2; fi
vmtouch -q -e "${files[@]}"
out=$(vmtouch -v "${files[@]}" 2>/dev/null | grep "Resident Pages" | tail -1)
pct=$(echo "$out" | sed -E 's/.*Resident Pages: *([0-9]+)\/([0-9]+).*/\1 \2/' | awk '{ if ($2 == 0) print 0; else printf "%.3f", 100*$1/$2 }')
echo "resident after eviction: ${pct}% ($out)"
awk -v p="$pct" 'BEGIN { exit (p < 1.0) ? 0 : 2 }'
