#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --time=36:00:00



ml load releases/2018b
ml load Python/3.6.6-foss-2018b
source $CECIHOME/venvsimGW/bin/activate


nrays=100
ri=1
rs=10
python create_linspace.py ray.in $ri $rs $nrays

num=100
bi=2
bs=8
python create_linspace.py mass.in $bi $bs $num

num2=1
di=25
ds=25
mkdir results
count2=1
while [ $count2 -le $nrays ]
do
    echo $count2
    paramIN2=$(cat ray.in | head -$count2| tail -1)
    srun -n1 --exclusive python carte_mass_sim.py $di $ds $num2 $paramIN2 $count2 &
    count2=$(($count2+1))
done
wait
deactivate
cd results
cat simRay* > simtotal.txt
