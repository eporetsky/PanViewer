#!/usr/bin/env bash
# Shared env bootstrap for Slurm/batch jobs.
# Sources conda + genetribe-upstream onto PATH. Fail loudly if genetribe missing.
#
# Optional overrides:
#   GENETRIBE_CONDA_ENV=/path/to/conda/envs/genetribe
#   GENETRIBE_HOME=/path/to/genetribe-upstream
#   CONDA_ROOT=/path/to/miniconda|mambaforge|anaconda
#
# shellcheck shell=bash

_genetribe_bootstrap() {
  local repo_root="${REPO_ROOT:?REPO_ROOT must be set}"
  local cand

  # --- conda env (blast, bedtools, jcvi) ---
  if [[ -z "${GENETRIBE_CONDA_ENV:-}" ]]; then
    for cand in \
      "${CONDA_PREFIX:-}" \
      "${HOME}/miniconda3/envs/genetribe" \
      "${HOME}/mambaforge/envs/genetribe" \
      "${HOME}/conda/envs/genetribe" \
      "${HOME}/.conda/envs/genetribe"
    do
      if [[ -n "${cand}" && -d "${cand}/bin" ]]; then
        GENETRIBE_CONDA_ENV="${cand}"
        break
      fi
    done
  fi

  if [[ -n "${GENETRIBE_CONDA_ENV:-}" && -d "${GENETRIBE_CONDA_ENV}/bin" ]]; then
    export PATH="${GENETRIBE_CONDA_ENV}/bin:${PATH}"
  fi

  # Prefer `conda activate` when possible so activate.d hooks run
  if [[ -z "${CONDA_ROOT:-}" ]]; then
    for cand in \
      "${HOME}/miniconda3" \
      "${HOME}/mambaforge" \
      "${HOME}/miniforge3" \
      "${HOME}/anaconda3" \
      "${HOME}/conda"
    do
      if [[ -f "${cand}/etc/profile.d/conda.sh" ]]; then
        CONDA_ROOT="${cand}"
        break
      fi
    done
  fi
  if [[ -n "${CONDA_ROOT:-}" && -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${CONDA_ROOT}/etc/profile.d/conda.sh"
    if [[ -n "${GENETRIBE_CONDA_ENV:-}" ]]; then
      conda activate "${GENETRIBE_CONDA_ENV}" 2>/dev/null \
        || conda activate genetribe 2>/dev/null \
        || true
    else
      conda activate genetribe 2>/dev/null || true
    fi
  fi

  # --- GeneTribe CLI (NOT provided by conda; comes from make setup clone) ---
  if [[ -z "${GENETRIBE_HOME:-}" ]]; then
    for cand in \
      "${repo_root}/genetribe-upstream" \
      "${GENETRIBE_CONDA_ENV:-}/../genetribe-upstream" \
      "${repo_root}"
    do
      if [[ -e "${cand}/genetribe" ]]; then
        GENETRIBE_HOME="${cand}"
        break
      fi
    done
  fi

  if [[ -n "${GENETRIBE_HOME:-}" && -e "${GENETRIBE_HOME}/genetribe" ]]; then
    export GENETRIBE_HOME
    export PATH="${GENETRIBE_HOME}:${PATH}"
  fi

  if ! command -v genetribe >/dev/null 2>&1; then
    echo "ERROR: genetribe not found on PATH." >&2
    echo "  Conda supplies blast/bedtools/jcvi only." >&2
    echo "  From this repository root run:  make setup" >&2
    echo "  That clones genetribe-upstream/ and runs install.sh." >&2
    echo "  REPO_ROOT=${repo_root}" >&2
    echo "  GENETRIBE_HOME=${GENETRIBE_HOME:-unset}" >&2
    echo "  GENETRIBE_CONDA_ENV=${GENETRIBE_CONDA_ENV:-unset}" >&2
    echo "  PATH=${PATH}" >&2
    return 1
  fi

  echo "bootstrap: genetribe=$(command -v genetribe)"
  echo "bootstrap: GENETRIBE_HOME=${GENETRIBE_HOME:-}"
  echo "bootstrap: GENETRIBE_CONDA_ENV=${GENETRIBE_CONDA_ENV:-}"
  echo "bootstrap: blast=$(command -v blastp || echo MISSING)"
  echo "bootstrap: bedtools=$(command -v bedtools || echo MISSING)"
  return 0
}

# Optional Slurm site flags for submit scripts.
# Set GT_ACCOUNT / GT_PARTITION if your cluster requires them.
_genetribe_sbatch_site_args() {
  SBATCH_SITE_ARGS=()
  if [[ -n "${GT_ACCOUNT:-}" ]]; then
    SBATCH_SITE_ARGS+=(--account="${GT_ACCOUNT}")
  fi
  if [[ -n "${GT_PARTITION:-}" ]]; then
    SBATCH_SITE_ARGS+=(--partition="${GT_PARTITION}")
  fi
}
