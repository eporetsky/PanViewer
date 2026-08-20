#!/usr/bin/env bash
# Put conda + GeneTribe CLI on PATH. Fail loudly if `genetribe` is missing.
#
# Requires REPO_ROOT to be set by the caller.
# Prefer: conda activate genetribe  (after `make setup`)
# Optional overrides: GENETRIBE_HOME, GENETRIBE_CONDA_ENV, CONDA_ROOT
# shellcheck shell=bash

_genetribe_bootstrap() {
  local repo_root="${REPO_ROOT:?REPO_ROOT must be set}"
  local cand

  # --- conda env (blast, bedtools, jcvi from environment.yml) ---
  if [[ -z "${GENETRIBE_CONDA_ENV:-}" ]]; then
    for cand in \
      "${CONDA_PREFIX:-}" \
      "${HOME}/conda/envs/genetribe" \
      "${HOME}/miniconda3/envs/genetribe" \
      "${HOME}/mambaforge/envs/genetribe"
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

  if [[ -z "${CONDA_ROOT:-}" ]]; then
    for cand in \
      "${HOME}/miniconda3" \
      "${HOME}/mambaforge" \
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

  # --- GeneTribe CLI (from make setup → genetribe-upstream/; not a conda package) ---
  if [[ -z "${GENETRIBE_HOME:-}" ]]; then
    for cand in \
      "${repo_root}/genetribe-upstream" \
      "${GENETRIBE_CONDA_ENV:-}/../genetribe-upstream"
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
    echo "  From genetribe/ run:  make setup && conda activate genetribe" >&2
    echo "  (environment.yml installs blast/bedtools/jcvi; make setup clones the CLI.)" >&2
    echo "  REPO_ROOT=${repo_root}" >&2
    echo "  GENETRIBE_HOME=${GENETRIBE_HOME:-unset}" >&2
    return 1
  fi

  echo "bootstrap: genetribe=$(command -v genetribe)"
  echo "bootstrap: blast=$(command -v blastp || echo MISSING)"
  echo "bootstrap: bedtools=$(command -v bedtools || echo MISSING)"
  return 0
}
