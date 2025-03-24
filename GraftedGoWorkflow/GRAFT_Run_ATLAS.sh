#!/bin/bash
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
NRESIDUE=$(awk '$1=="ATOM"&&$3=="CA"{print $6}' $1 | tail -n 1)

awk -v nres=$NRESIDUE 'NR>1&&$1<nres{m=m+1}END{print m}' mask.dat >  $simul/mask.dat
awk -v nres=$NRESIDUE 'NR>1&&$1<nres{print $1}' mask.dat >> $simul/mask.dat

awk -v nres=$NRESIDUE 'NR>1&&$1<nres{m=m+1}END{print m}' mask2.dat >  $simul/mask2.dat
awk -v nres=$NRESIDUE 'NR>1&&$1<nres{print $1}' mask2.dat >> $simul/mask2.dat


#cp ignore.dat $simul
cp ATLAS.sh $simul
echo $seed
cd $simul
place=$(pwd)
echo $place

##################### SET PARAMS #############################
sed -i "s+Seed=.*#+Seed= $seed #+g" param_input.dat
sed -i "s+Temp=.*#+Temp= $tempera #+g" param_input.dat
sed -i "s+cd FOLDER+cd $place+g" ATLAS.sh
sed -i "s+INFILE+in.${tipo}_lmp+g" ATLAS.sh
sed -i "s+NAME+${tipo}-P-${seed}-T-${tempera}+g" ATLAS.sh
################### CREATE SIMULATION FILES ##################
contactmap-Shae_NY_design-GO_GRAFTFLOW-COMPLEX.csh $1 param_input.dat
GO_PROTGRAFT.bash $tipo"-out"
################### SUBMIT JOB #######################
sbatch ATLAS.sh > jobid.dat
cd $start
done
done
