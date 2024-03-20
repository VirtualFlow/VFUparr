#!/bin/bash
#SBATCH --account=TODO # Please provde the username
#SBATCH --ntasks-per-node=20
#SBATCH --mem=7000M               # memory (per node)
#SBATCH --time=12:0:00
#SBATCH --job-name VFUparr
#SBATCH --array=1-999
#SBATCH -e stderr.txt
#SBATCH -o stdout.txt
#SBATCH --open-mode=append
#SBATCH --export=NONE

module --force purge
module load nixpkgs/16.09
module load gcc/7.3.0
module load rdkit/2019.03.4
module load scipy-stack/2019b
module load openbabel

python3 dataset_calc.py $SLURM_ARRAY_TASK_ID
