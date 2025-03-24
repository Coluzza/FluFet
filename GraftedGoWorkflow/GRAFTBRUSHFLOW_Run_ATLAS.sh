#!/bin/bash
#FIRST ARGUMENT PDB_FILE SECOND ARGUMENT PARAM_FILE
tipo=$(echo $1| awk -F ".pdb" '{print $1}')
latest=$(echo $tipo/Simul_`date +%d-%m_%H%M%S`)
start=$PWD
for seed in {1..5}
do
for tempera in {2..10..2}
do
for fraction in 0.3 0.5 0.9
do
simul=$latest"/P-"$seed"-T-"$tempera"-F-"$fraction
echo $simul

mkdir -p $simul
cp $1 $simul

cp $2 $simul/param_input.dat
NRESIDUE=$(awk '$1=="ATOM"&&$3=="CA"{print $6}' $1 | tail -n 1)

awk -v nres=$NRESIDUE 'NR>1&&$1<nres{m=m+1}END{print m}' mask.dat >  $simul/mask.dat
awk -v nres=$NRESIDUE 'NR>1&&$1<nres{print $1}' mask.dat >> $simul/mask.dat

awk -v nres=$NRESIDUE 'NR>1&&$1<nres{m=m+1}END{print m}' mask2.dat >  $simul/mask2.dat
awk -v nres=$NRESIDUE 'NR>1&&$1<nres{print $1}' mask2.dat >> $simul/mask2.dat


#cp ignore.dat $simul
cp ATLAS_FLOW.sh $simul
echo $seed
cd $simul
place=$(pwd)
echo $place

##################### SET PARAMS #############################
sed -i "s+Seed=.*#+Seed= $seed #+g" param_input.dat
sed -i "s+Temp=.*#+Temp= $tempera #+g" param_input.dat
sed -i "s+Polymer_Fraction=.*#+Polymer_Fraction= $fraction #+g" param_input.dat
sed -i "s+cd FOLDER+cd $place+g" ATLAS_FLOW.sh
sed -i "s+INFILE+in.${tipo}_lmp+g" ATLAS_FLOW.sh
sed -i "s+NAME+${tipo}-P-${seed}-T-${tempera}-T-${fraction}+g" ATLAS_FLOW.sh
################### CREATE SIMULATION FILES ##################
GRAFTFLOW-BRUSH-COMPLEX.sh $1 param_input.dat
GO_PROTGRAFT_BRUSH_FLOW.bash $tipo"-out"
################### SUBMIT JOB #######################
sbatch ATLAS_FLOW.sh > jobid.dat
cd $start
done
done
done
