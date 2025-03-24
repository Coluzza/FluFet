#!/bin/bash

tipo=$(echo $1| awk -F ".pdb" '{print $1}')
latest=$(echo $tipo/Simul_`date +%d-%m_%H%M%S`)
start=$PWD
for seed in {1..10}
do
simul=$latest"/P-"$seed
mkdir -p $simul
cp $1 $simul


cp $2 $simul/param_input.dat
cp mask.dat $simul
#cp ignore.dat $simul
cp GT4.sh $simul
echo $seed
cd $simul
place=$(pwd)
echo $place

##################### SET PARAMS #############################
sed -i "s+Seed= .* #+Seed= $seed #+g" param_input.dat
sed -i "s+cd FOLDER+cd $place+g" GT4.sh
sed -i "s+in.Rudimero_lmp+in.${tipo}_lmp+g" GT4.sh
################### CREATE SIMULATION FILES ##################
contactmap-Shae_NY_design-GO_GRAFTFLOW-COMPLEX.csh $1 param_input.dat
GO_PROTFLOW.bash $tipo"-out"
################### SUBMIT JOB #######################
sbatch GT4.sh > jobid.dat
cd $start
done

