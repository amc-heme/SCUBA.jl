# Regenerate derived fixtures and R goldens. Run from the repository root: make -C fixtures
DATA ?= data
PY   ?= conda run -n scuba-fixtures python

.PHONY: all derived goldens pack
all: derived goldens

derived: $(DATA)/AML_h5ad.h5ad
	cd .. && $(PY) fixtures/make_fixtures.py fixtures/$(DATA)

goldens:
	cd .. && Rscript fixtures/make_goldens.R fixtures/build/goldens

pack:
	cd .. && julia --project=. fixtures/pack.jl
