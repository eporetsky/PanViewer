SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c

.DEFAULT_GOAL := setup

GENETRIBE_REPO ?= https://github.com/chenym1/genetribe.git
GENETRIBE_DIR ?= genetribe-upstream
CONDA_ENV ?= genetribe
ENV_YML ?= environment.yml

.PHONY: setup path

setup:
	@command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found on PATH. Install Miniconda/Mambaforge and retry." >&2; exit 1; }
	@echo "Preparing local GeneTribe checkout..."
	@if [[ -d "$(GENETRIBE_DIR)/.git" ]]; then \
		echo "Updating $(GENETRIBE_DIR)..."; \
		git -C "$(GENETRIBE_DIR)" pull --ff-only; \
	elif [[ -e "$(GENETRIBE_DIR)" ]]; then \
		echo "ERROR: $(GENETRIBE_DIR) exists but is not a git clone. Remove or rename it, then rerun." >&2; \
		exit 1; \
	else \
		echo "Cloning $(GENETRIBE_REPO) -> $(GENETRIBE_DIR)..."; \
		git clone "$(GENETRIBE_REPO)" "$(GENETRIBE_DIR)"; \
	fi
	@echo "Running GeneTribe install.sh..."
	@cd "$(GENETRIBE_DIR)" && bash ./install.sh
	@test -x "$(GENETRIBE_DIR)/genetribe" || test -f "$(GENETRIBE_DIR)/genetribe" || { \
		echo "ERROR: $(GENETRIBE_DIR)/genetribe missing after install.sh." >&2; exit 1; }
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
		'export GENETRIBE_HOME="$(CURDIR)/$(GENETRIBE_DIR)"' \
		'export _GENETRIBE_OLD_PATH="$$PATH"' \
		'export PATH="$$GENETRIBE_HOME:$$PATH"' \
		> "$$ENV_PREFIX/etc/conda/activate.d/genetribe.sh"; \
	printf '%s\n' \
		'if [[ -n "$${_GENETRIBE_OLD_PATH:-}" ]]; then' \
		'  export PATH="$$_GENETRIBE_OLD_PATH"' \
		'  unset _GENETRIBE_OLD_PATH' \
		'fi' \
		'unset GENETRIBE_HOME' \
		> "$$ENV_PREFIX/etc/conda/deactivate.d/genetribe.sh"
	@echo ""
	@echo "Setup complete."
	@echo "Run: conda activate $(CONDA_ENV)"
	@echo ""
	@printf '%s\n' "Note: genetribe lives in ./$(GENETRIBE_DIR) (added to PATH via activate.d)."
	@printf '%s\n' "If 'genetribe: command not found' after conda activate:"
	@printf '%s\n' "  conda deactivate && conda activate $(CONDA_ENV)"
	@printf '%s\n' "Or paste from this directory:"
	@printf '  export PATH="%s:$$PATH"\n\n' "$(CURDIR)/$(GENETRIBE_DIR)"

path:
	@printf '%s\n' "From genetribe/ directory:"
	@printf '  export PATH="%s:$$PATH"\n' "$(CURDIR)/$(GENETRIBE_DIR)"
