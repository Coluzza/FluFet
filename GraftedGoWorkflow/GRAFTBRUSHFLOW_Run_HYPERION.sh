#!/bin/bash
Prepare_HYPERION () {
		
cat << EOF > HYPERION.sh

#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --partition=regular
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:rtx3090:1

source ~/.bashrc

module purge

module load LAMMPS/23Jun2022-foss-2021b-kokkos-CUDA-11.4.1
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.1pgb-out >> OUT 2>&1

}

Prepare_HYPERION_FLOW () {
		
cat << EOF > HYPERION_FLOW.sh

#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$3
#SBATCH --partition=regular
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=10
#SBATCH --gres=gpu:rtx3090:1
#SBATCH -d afterany:$2
source ~/.bashrc

module purge

module load LAMMPS/23Jun2022-foss-2021b-kokkos-CUDA-11.4.1
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 10 lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.1pgb-srd  >> OUT 2>&1

}


#ARGUMENTS     #1       #2        #3
#          PDB_FILE PARAM_FILE PKA_FILE

pdb=realpath $1
param=realpath $2
pka=realpath $3

proteina=$(echo $pdb|xargs -n 1 basename | awk -F ".pdb" '{print $1}')
geometry=$(grep Geometry $param | awk '{print $2}')
type=$(grep Simul_Type $param | awk '{print $2}')

NAME=$(echo $proteina-$geometry-$type)



mkdir -p $type-$geometry/$proteina
cd $type-$geometry/$proteina
start=$PWD

latest=$(echo Simul_`date +%d-%m_%H%M%S`)

for seed in {1..5}
do
for tempera in {2..10..2}
do
for fraction in 0.3 0.5 0.9
do
for ph in 7.0
do

simul=$latest"/P-"$seed"-T-"$tempera"-F-"$fraction"-pH-"$ph
mkdir -p $simul

cd $simul
place=$(pwd)
echo $place

cp $pdb ./

cp $param ./param_input.dat

##################### SET PARAMS #############################
sed -i "s+Seed=.*#+Seed= $seed #+g" param_input.dat
sed -i "s+Temp=.*#+Temp= $tempera #+g" param_input.dat
sed -i "s+Polymer_Fraction=.*#+Polymer_Fraction= $fraction #+g" param_input.dat

################### CREATE SIMULATION FILES ##################
GENERAL_PREP.sh $pdb param_input.dat $pka $ph
################### SUBMIT JOB #######################
Prepare_HYPERION($pwd,$NAME)
sbatch HYPERION.sh > jobid.dat
job=$(echo jobid.dat)
Prepare_HYPERION_FLOW($pwd,$job,$NAME)
sbatch HYPERION_FLOW.sh > jobid_flow.dat
cd $start
done
done
done
done




