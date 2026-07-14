SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c

.DEFAULT_GOAL := setup

PANDAGMA_REPO ?= https://github.com/legumeinfo/pandagma.git
PANDAGMA_DIR ?= pandagma-upstream
CONDA_ENV ?= pandagma
ENV_YML ?= environment.yml

.PHONY: setup path patch-pandagma-ingest-fasta patch-pandagma-dagchainer

# Ensure ingest writes 02_fasta_nuc/<stem>.$fa (upstream omits suffix → broken globs / awk).
patch-pandagma-ingest-fasta:
	@bash "$(CURDIR)/scripts/patch-pandagma-pan-ingest-fasta.sh"

patch-pandagma-dagchainer:
	@bash "$(CURDIR)/scripts/patch-pandagma-pan-dagchainer.sh"

setup:
	@command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found on PATH. Install Miniconda/Mambaforge and retry." >&2; exit 1; }
	@echo "Preparing local pandagma checkout..."
	@if [[ -d "$(PANDAGMA_DIR)/.git" ]]; then \
		echo "Updating $(PANDAGMA_DIR)..."; \
		git -C "$(PANDAGMA_DIR)" pull --ff-only; \
	elif [[ -e "$(PANDAGMA_DIR)" ]]; then \
		echo "ERROR: $(PANDAGMA_DIR) exists but is not a git clone. Remove or rename it, then rerun." >&2; \
		exit 1; \
	else \
		echo "Cloning $(PANDAGMA_REPO) -> $(PANDAGMA_DIR)..."; \
		git clone "$(PANDAGMA_REPO)" "$(PANDAGMA_DIR)"; \
	fi
	@test -x "$(PANDAGMA_DIR)/bin/pandagma" || { echo "ERROR: $(PANDAGMA_DIR)/bin/pandagma missing or not executable." >&2; exit 1; }
	@echo "Patching $(PANDAGMA_DIR)/bin/pandagma-pan.sh (ingest FASTA suffix fix)..."
	@bash "$(CURDIR)/scripts/patch-pandagma-pan-ingest-fasta.sh"
	@echo "Patching $(PANDAGMA_DIR)/bin/pandagma-pan.sh (dagchainer batch reliability)..."
	@bash "$(CURDIR)/scripts/patch-pandagma-pan-dagchainer.sh"
	@echo "Creating/updating conda env $(CONDA_ENV)..."
	@if conda env list | grep -qE '^$(CONDA_ENV)[[:space:]]'; then \
		conda env update -f "$(ENV_YML)" --prune -n "$(CONDA_ENV)"; \
	else \
		conda env create -f "$(ENV_YML)"; \
	fi
	@ENV_PREFIX="$$(conda run -n "$(CONDA_ENV)" bash -c 'echo "$$CONDA_PREFIX"')"; \
	test -n "$$ENV_PREFIX" || { echo "ERROR: could not resolve conda prefix for env $(CONDA_ENV) (conda run failed?)" >&2; exit 1; }; \
	mkdir -p "$$ENV_PREFIX/etc/conda/activate.d" "$$ENV_PREFIX/etc/conda/deactivate.d"; \
	printf '%s\n' \
		'export PANDAGMA_HOME="$(CURDIR)/$(PANDAGMA_DIR)"' \
		'export PANDAGMA_BIN="$$PANDAGMA_HOME/bin"' \
		'export _PANDAGMA_OLD_PATH="$$PATH"' \
		'export PATH="$$PANDAGMA_BIN:$$PATH"' \
		> "$$ENV_PREFIX/etc/conda/activate.d/pandagma.sh"; \
	printf '%s\n' \
		'if [[ -n "$${_PANDAGMA_OLD_PATH:-}" ]]; then' \
		'  export PATH="$$_PANDAGMA_OLD_PATH"' \
		'  unset _PANDAGMA_OLD_PATH' \
		'fi' \
		'unset PANDAGMA_HOME' \
		'unset PANDAGMA_BIN' \
		> "$$ENV_PREFIX/etc/conda/deactivate.d/pandagma.sh"
	@echo ""
	@echo "Setup complete."
	@echo "Run: conda activate $(CONDA_ENV)"
	@echo ""
	@printf '%s\n' "Note: pandagma is not a conda-installed binary; PATH must include ./$(PANDAGMA_DIR)/bin."
	@printf '%s\n' "If 'pandagma: command not found' after conda activate:"
	@printf '%s\n' "  conda deactivate && conda activate $(CONDA_ENV)"
	@printf '%s\n' "Or paste from repo root:"
	@printf '  export PATH="%s/bin:$$PATH"\n\n' "$(CURDIR)/$(PANDAGMA_DIR)"

path:
	@printf '%s\n' "From project root:"
	@printf '  export PATH="%s/bin:$$PATH"\n' "$(CURDIR)/$(PANDAGMA_DIR)"
