.PHONY: clean clean_all readme benchmarks correctness

PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

EXTENSION_NAME=duckvep

# Set to 1 to enable Unstable API (binaries will only work on TARGET_DUCKDB_VERSION, forwards compatibility will be broken)
# Note: currently extension-template-rs requires this, as duckdb-rs relies on unstable C API functionality
USE_UNSTABLE_C_API=1

# Target DuckDB version
TARGET_DUCKDB_VERSION=v1.5.3

all: configure debug

# Include makefiles from DuckDB
include extension-ci-tools/makefiles/c_api_extensions/base.Makefile
include extension-ci-tools/makefiles/c_api_extensions/rust.Makefile

configure: venv platform extension_version

debug: build_extension_library_debug build_extension_with_metadata_debug
release: build_extension_library_release build_extension_with_metadata_release

test: test_debug
test_debug: test_extension_debug
test_release: test_extension_release

clean: clean_build clean_rust
clean_all: clean_configure clean

# Rendered docs are generated from their .Rmd + the data files they import, so the output is a
# real file target with explicit prerequisites: editing the .Rmd, the render script, or any
# imported CSV invalidates the rendered .md and `make` re-runs only what is stale. The phony
# names (readme/benchmarks/correctness) remain as convenience aliases. README.md also embeds
# live SQL via duckknit, so it additionally needs the built extension (`make debug`).

README.md: README.Rmd scripts/render-readme.R \
           build/debug/duckvep.duckdb_extension \
           correctness/data/concordance_by_impact.csv \
           correctness/data/error_transitions.csv \
           benchmarks/data/timings.csv
	Rscript scripts/render-readme.R
readme: README.md

# render-correctness.R renders both correctness.md and cache-build/README.md.
correctness/correctness.md: correctness/correctness.Rmd correctness/cache-build/README.Rmd \
           scripts/render-correctness.R \
           correctness/data/concordance_by_impact.csv \
           correctness/data/discordance_by_consequence.csv \
           correctness/data/error_transitions.csv \
           correctness/data/methodology_audit.csv \
           correctness/data/so_term_transitions.csv \
           correctness/cache-build/data/cache_stats.csv
	Rscript scripts/render-correctness.R
correctness: correctness/correctness.md

benchmarks/results.md: benchmarks/results.Rmd scripts/render-benchmarks.R \
           benchmarks/data/footprint.csv \
           benchmarks/data/thread_scaling.csv \
           benchmarks/data/throughput.csv \
           benchmarks/data/timings.csv
	Rscript scripts/render-benchmarks.R
benchmarks: benchmarks/results.md
