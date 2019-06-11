#!/bin/bash
#SBATCH -J LHSSim
#SBATCH --cpus-per-task=6
#SBATCH --partition=med
#SBATCH --mem=5Gb
#SBATCH --time=03:00:00
#SBATCH --mail-user=mwilliam@ucdavis.edu
#SBATCH --mail-type=ALL

SCRATCH=/scratch/$USER/$SLURM_JOB_ID
echo $SLURM_SUBMIT_DIR
echo Creating temp dir $SCRATCH
srun mkdir -p $SCRATCH || exit $?
srun cp -r $SLURM_SUBMIT_DIR/Scripts/  $SCRATCH || exit $?

module load R
hostname -f
echo $SLURM_ARRAY_TASK_ID

cd $SCRATCH
srun ls
Rscript --vanilla ./Scripts/onerun.R ${1} ${2}

srun cp -r $SCRATCH/*.rds $SLURM_SUBMIT_DIR/Outs/

echo Removing $SCRATCH
srun rm -rf $SCRATCH || exit $?
