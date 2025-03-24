#!/bin/bash
NUMBER_PROC=10

#FIRST ARGUMENT PDB_FILE SECOND ARGUMENT PARAM_FILE
tipo=$(echo $1| awk -F ".pdb" '{print $1}')
latest=$(echo $tipo/Simul_`date +%d-%m_%H%M%S`)
start=$PWD
for seed in {1..10}
do
for tempera in {2..10}
do
simul=$latest"/P-"$seed"-T-"$tempera
mkdir -p $simul
cp $1 $simul


cp $2 $simul/param_input.dat
cp mask.dat $simul
cp mask2.dat $simul
#cp ignore.dat $simul
cp ATLAS.sh $simul
echo $seed
cd $simul
place=$(pwd)
echo $place

##################### SET PARAMS #############################
sed -i "s+Seed= .* #+Seed= $seed #+g" param_input.dat
sed -i "s+Temp= .* #+Temp= $tempera #+g" param_input.dat
sed -i "s+cd FOLDER+cd $place+g" ATLAS.sh
sed -i "s+INFILE+in.${tipo}_lmp+g" ATLAS.sh
################### CREATE SIMULATION FILES ##################
contactmap-Shae_NY_design-GO_GRAFTFLOW-COMPLEX.csh $1 param_input.dat
GO_PROTGRAFT.bash $tipo"-out"
################### SUBMIT JOB #######################
mpirun -np $NUMBER_PROC lmp_mpi -i in.${tipo}_lmp >> OUT 2>&1

cd $start
done
done
