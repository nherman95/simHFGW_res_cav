#!/bin/bash

ml releases/2018b
ml Python/3.6.6-foss-2018b
ml GCC
source $CECIHOME/venvsimGW/bin/activate
make
num=50
bi=2
bs=8
python create_linspace.py mass.in $bi $bs $num

num2=50
di=15
ds=30
mkdir results
cat <<EOM >"jobarray.sh"
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
source $CECIHOME/venvsimGW/bin/activate
paramIN=$(echo '$(cat mass.in | head -$SLURM_ARRAY_TASK_ID| tail -1)')
python mass_sim.py $(echo '$SLURM_ARRAY_TASK_ID $paramIN') $di $ds $num2
deactivate
EOM
sbatch jobarray.sh
