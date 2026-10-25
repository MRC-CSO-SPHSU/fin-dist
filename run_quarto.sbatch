#!/bin/bash
#SBATCH --job-name=fd_targets
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/targets_%j.out
#SBATCH --error=logs/targets_%j.err
#SBATCH --account=none
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=darwin.delcastillofernandez@glasgow.ac.uk

set -uo pipefail

module purge
module load apps/miniforge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate quarto

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export RANGER_NUM_THREADS=1

# Node-local scratch for worker spillover and targets temp files.
export TMPDIR="${SLURM_TMPDIR:-/tmp}"

echo "Job ${SLURM_JOB_ID} on $(hostname)"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"
echo "R: $(which R)"; R --version | head -n 1
echo "Quarto: $(which quarto)"; quarto --version

cd "${SLURM_SUBMIT_DIR}"

# tar_make_future() dispatches each worker target as its own SLURM job via future.batchtools
# The final tar_quarto(report) target renders 06_full_models.qmd.
Rscript -e 'targets::tar_make_future(workers = 16L)' 2>&1 | tee "targets_make_${SLURM_JOB_ID}.log"
status=${PIPESTATUS[0]}
echo "tar_make_future exit status: $status"

exit $status
