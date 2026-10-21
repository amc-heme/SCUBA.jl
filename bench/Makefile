# SCUBA benchmark suite. Run from this directory. Datasets live in $(DATA) (default ~/scuba-bench-data).
#   make preflight            check tools, disk, RAM; write results/env/<RUN_ID>.json
#   make smoke                10k end-to-end (generate, Julia + Python harness, merge, report); R=1 adds the R harness
#   make all-100k / all-1M    full pass at one scale (data, R objects, all harnesses, cold, concurrency, merge, report)
#   make data-3M              opt-in 3M dataset (preflight first)
SHELL := /bin/bash
CONDA_ENV ?= scuba-fixtures
CONDA := source ~/miniconda3/etc/profile.d/conda.sh && conda run -n $(CONDA_ENV) --no-capture-output
PY := $(CONDA) python -I
RUN_ID ?= $(shell date -u +%Y%m%dT%H%M%S)
DATA ?= $(HOME)/scuba-bench-data
JULIA_THREADS ?= 1
TOOLPY := python3 -I
TIMEOUT := $(TOOLPY) tools/run_with_timeout.py --timeout 14400 --
export JULIA_CONDAPKG_BACKEND=Null
export PYTHON_JULIAPKG_PROJECT=$(abspath ..)/python/.juliapkg-test

.PHONY: preflight smoke merge report clean-results
preflight:
	$(PY) tools/preflight.py --run-id $(RUN_ID) $(if $(SCALE),--scale $(SCALE),)

data-%:
	$(PY) tools/preflight.py --scale $* --run-id $(RUN_ID) --no-env
	$(PY) generate/make_dataset.py --scale $* --data-dir $(DATA)

r-objects-%:
	Rscript generate/make_r_objects.R --scale $* --data-dir $(DATA)

julia-%:
	$(TIMEOUT) julia --project=. -t $(JULIA_THREADS) julia/run.jl --scale $* --run-id $(RUN_ID) --resume $(JULIA_FLAGS)

r-%:
	$(TIMEOUT) Rscript r/run.R --scale $* --engine r --run-id $(RUN_ID) --resume $(R_FLAGS)
	$(TIMEOUT) Rscript r/run.R --scale $* --engine julia --run-id $(RUN_ID) --resume $(R_FLAGS)

python-%:
	$(TIMEOUT) $(PY) python/run.py --scale $* --run-id $(RUN_ID) --resume $(PY_FLAGS)

cold-%:
	for impl in jl_h5ad_csc jl_bpcells jl_bpcells_rowmajor; do $(TOOLPY) tools/run_cold.py --scale $* --lang julia --impl $$impl --run-id $(RUN_ID); done
	for impl in r_seurat_bpcells r_seurat_bpcells_rowmajor r_sce_backed; do $(TOOLPY) tools/run_cold.py --scale $* --lang r --impl $$impl --engine r --run-id $(RUN_ID); done
	$(TOOLPY) tools/run_cold.py --scale $* --lang r --impl r_seurat_bpcells --engine julia --run-id $(RUN_ID)
	$(TOOLPY) tools/run_cold.py --scale $* --lang python --impl raw_anndata_py --run-id $(RUN_ID)

conc-%:
	$(TOOLPY) tools/concurrency.py --scale $* --lang julia --impl jl_h5ad_csc --run-id $(RUN_ID)
	$(TOOLPY) tools/concurrency.py --scale $* --lang julia --impl jl_bpcells --run-id $(RUN_ID)
	$(TOOLPY) tools/concurrency.py --scale $* --lang julia --impl jl_h5ad_csc --mode threads --share shared --run-id $(RUN_ID)
	$(TOOLPY) tools/concurrency.py --scale $* --lang julia --impl jl_bpcells --mode threads --share per_task --run-id $(RUN_ID)
	$(TOOLPY) tools/concurrency.py --scale $* --lang julia --impl jl_bpcells_rowmajor --run-id $(RUN_ID)
	$(TOOLPY) tools/concurrency.py --scale $* --lang r --impl r_seurat_bpcells_rowmajor --engine r --run-id $(RUN_ID) --sessions 1,2,4,8
	$(TOOLPY) tools/concurrency.py --scale $* --lang r --impl r_seurat_bpcells_rowmajor --engine julia --run-id $(RUN_ID) --sessions 1,2,4,8
	$(TOOLPY) tools/concurrency.py --scale $* --lang r --impl r_seurat_bpcells --engine r --run-id $(RUN_ID) --sessions 1,2,4,8
	$(TOOLPY) tools/concurrency.py --scale $* --lang r --impl r_seurat_bpcells --engine julia --run-id $(RUN_ID) --sessions 1,2,4,8
	$(TOOLPY) tools/concurrency.py --scale $* --lang python --impl raw_anndata_py --run-id $(RUN_ID) --sessions 1,2,4,8

merge:
	$(PY) tools/merge_results.py

report: merge
	Rscript report/make_report.R
	-Rscript -e 'rmarkdown::render("report/report.Rmd", quiet = TRUE)'

all-%: data-% r-objects-% julia-% python-% r-% cold-% conc-% report
	@echo "done: $* (run $(RUN_ID))"

smoke: RUN_ID := smoke
smoke:
	$(PY) tools/preflight.py --scale 10k --run-id smoke
	$(PY) generate/make_dataset.py --scale 10k --data-dir $(DATA)
	$(TIMEOUT) julia --project=. -t 1 julia/run.jl --scale 10k --run-id smoke --only jl_h5ad_csc,jl_h5ad_csr
	$(TIMEOUT) $(PY) python/run.py --scale 10k --run-id smoke --only raw_anndata_py,raw_h5py_py
ifdef R
	Rscript generate/make_r_objects.R --scale 10k --data-dir $(DATA)
	$(TIMEOUT) julia --project=. -t 1 julia/run.jl --scale 10k --run-id smoke --only jl_bpcells
	$(TIMEOUT) Rscript r/run.R --scale 10k --engine r --run-id smoke
	$(TIMEOUT) Rscript r/run.R --scale 10k --engine julia --run-id smoke
endif
	$(PY) tools/merge_results.py --run-id smoke
	Rscript report/make_report.R

clean-results:
	rm -rf results/raw results/results.parquet results/results.csv report/figures report/summary.md report/report.html
