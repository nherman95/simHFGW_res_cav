#!/bin/bash

ml releases/2018b
ml Python/3.6.6-foss-2018b
source $CECIHOME/venvsimGW/bin/activate


nrays=400
ri=1
rs=10
python create_linspace.py ray.in $ri $rs $nrays

num=400
bi=2
bs=8
python create_linspace.py mass.in $bi $bs $num

num2=1
di=25
ds=25
mkdir results
cat <<EOM >"jobarray.sh"
#!/bin/bash
#SBATCH --job-name=simmmass
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=1024
#SBATCH --array=1-400

ml load releases/2018b
ml load Python/3.6.6-foss-2018b
source $CECIHOME/venvsimGW/bin/activate
paramIN2=$(echo '$(cat ray.in | head -$SLURM_ARRAY_TASK_ID| tail -1)')
echo 'python bodemass_par.py $di $ds $num2 $paramIN2 $SLURM_ARRAY_TASK_ID'
deactivate
EOM
sbatch jobarray.sh
