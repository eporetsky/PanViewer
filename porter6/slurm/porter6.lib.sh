#!/usr/bin/env bash
# Shared by porter6.sh (sbatch) and porter6.batch.sh (login-node submitter).
# Run one FASTA through the pipeline; results land in results/<stem>/.
#
# Parallel jobs: set PORTER6_ISOLATE=1 so each job rsyncs a private copy of the
# pipeline (see PORTER6_TEMPLATE / PORTER6_STAGING). Otherwise jobs share porter6/ and collide.

_CPUS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="${_CPUS}"
export MKL_NUM_THREADS="${_CPUS}"
export OPENBLAS_NUM_THREADS="${_CPUS}"
export NUMEXPR_NUM_THREADS="${_CPUS}"
unset _CPUS

porter6_stem_from_fasta() {
  local _base
  _base=$(basename "$1")
  local STEM="$_base"
  local _ext
  for _ext in .fasta .FASTA .faa .FAA .fa .FA; do
    if [[ "$STEM" == *"$_ext" ]]; then
      STEM="${STEM%"$_ext"}"
      break
    fi
  done
  printf '%s' "$STEM"
}

porter6_run_one() {
  local FASTA="$1"
  local WORKDIR="${WORKDIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
  cd "$WORKDIR"

  if [[ "$FASTA" != /* ]]; then
    FASTA="${WORKDIR}/${FASTA}"
  fi
  if [[ ! -f "$FASTA" ]]; then
    echo "Error: FASTA not found: $FASTA" >&2
    return 1
  fi

  local STEM
  STEM=$(porter6_stem_from_fasta "$FASTA")

  local PORTER6_ROOT="${PORTER6_ROOT:-${WORKDIR}/porter6}"
  if [[ "${PORTER6_ISOLATE:-0}" == "1" ]]; then
    local TEMPLATE="${PORTER6_TEMPLATE:-${WORKDIR}/porter6}"
    local JOBTAG="${SLURM_JOB_ID:-$$}"
    local STAGING_ROOT="${PORTER6_STAGING:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}"
    local STAGING="${STAGING_ROOT}/porter6_work_${JOBTAG}"
    echo "PORTER6_ISOLATE: rsync ${TEMPLATE} -> ${STAGING}"
    rm -rf "${STAGING}"
    mkdir -p "${STAGING}"
    rsync -a "${TEMPLATE}/" "${STAGING}/"
    PORTER6_ROOT="${STAGING}"
  fi

  local RESULTS_DIR="${WORKDIR}/results/${STEM}"
  mkdir -p "$RESULTS_DIR" "${WORKDIR}/log"
  mkdir -p "${PORTER6_ROOT}/data/features/esm2"

  rm -f \
    "${PORTER6_ROOT}/predictor3/test_ensemble.json" \
    "${PORTER6_ROOT}/predictor3/output_predictions.csv" \
    "${PORTER6_ROOT}/predictor8/test_ensemble.json" \
    "${PORTER6_ROOT}/predictor8/output_predictions.csv" \
    "${PORTER6_ROOT}/porter6_predictions_merged.csv"
  find "${PORTER6_ROOT}/data/features/esm2" -maxdepth 1 -name '*.npy' -delete 2>/dev/null || true

  cp -f "$FASTA" "${PORTER6_ROOT}/data/dataset/set3.fasta"

  date

  # Site-specific: load CUDA if your cluster uses environment modules
  if [[ -n "${PORTER6_CUDA_MODULE:-}" ]]; then
    module load "${PORTER6_CUDA_MODULE}"
  fi

  export MPLBACKEND=Agg
  export TORCH_HOME="${TORCH_HOME:-$HOME/models/torch_models}"
  export HF_HOME="${HF_HOME:-$HOME/models/huggingface}"

  python -c "import torch; print('CUDA available:', torch.cuda.is_available()); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu')"

  chmod +x "${PORTER6_ROOT}/porter6.sh"
  "${PORTER6_ROOT}/porter6.sh"

  python "${PORTER6_ROOT}/merge_q3_q8.py"

  [[ -f "${PORTER6_ROOT}/predictor3/output_predictions.csv" ]] && mv -f "${PORTER6_ROOT}/predictor3/output_predictions.csv" "${RESULTS_DIR}/output_predictions_q3.csv"
  [[ -f "${PORTER6_ROOT}/predictor8/output_predictions.csv" ]] && mv -f "${PORTER6_ROOT}/predictor8/output_predictions.csv" "${RESULTS_DIR}/output_predictions_q8.csv"
  [[ -f "${PORTER6_ROOT}/predictor3/test_ensemble.json" ]] && mv -f "${PORTER6_ROOT}/predictor3/test_ensemble.json" "${RESULTS_DIR}/test_ensemble_q3.json"
  [[ -f "${PORTER6_ROOT}/predictor8/test_ensemble.json" ]] && mv -f "${PORTER6_ROOT}/predictor8/test_ensemble.json" "${RESULTS_DIR}/test_ensemble_q8.json"
  [[ -f "${PORTER6_ROOT}/porter6_predictions_merged.csv" ]] && mv -f "${PORTER6_ROOT}/porter6_predictions_merged.csv" "${RESULTS_DIR}/porter6_predictions_merged.csv"

  echo "Results in: ${RESULTS_DIR}"

  date
}
