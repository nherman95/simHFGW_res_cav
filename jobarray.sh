#!/bin/bash
#SBATCH --job-name=mass_sim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --array=1-50

ml load releases/2018b
ml load Python/3.6.6-foss-2018b
source /venvsimGW/bin/activate
paramIN=$(cat mass.in | head -$SLURM_ARRAY_TASK_ID| tail -1)
python mass_sim.py $SLURM_ARRAY_TASK_ID $paramIN 15 30 50
deactivate
