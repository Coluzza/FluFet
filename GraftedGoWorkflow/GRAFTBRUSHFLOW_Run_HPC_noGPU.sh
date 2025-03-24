#!/bin/bash
Prepare_HYPERION_LOWT () {
		
cat << EOF > HYPERION_LOWT.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --partition=regular
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=2


source ~/.bashrc

module purge

module load LAMMPS/23Jun2022-foss-2021b-kokkos-CUDA-11.4.1
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi_nogpu -i $1/in.$3-out-lowT >> OUT 2>&1
EOF
}

Prepare_HYPERION_HIGHT () {
		
cat << EOF > HYPERION_HIGHT.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --partition=regular
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=2

#SBATCH -d afterany:$4

source ~/.bashrc

module purge

module load LAMMPS/23Jun2022-foss-2021b-kokkos-CUDA-11.4.1
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi_nogpu -i $1/in.$3-out-highT >> OUT 2>&1
EOF
}

Prepare_HYPERION_FLOW () {
		
cat << EOF > HYPERION_FLOW.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --partition=regular
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=10

#SBATCH -d afterany:$4
source ~/.bashrc

module purge

module load LAMMPS/23Jun2022-foss-2021b-kokkos-CUDA-11.4.1
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 10 lmp_mpi_nogpu -i $1/in.$3-srd  >> OUT 2>&1
EOF
}
Prepare_NOTS_LOWT () {
cat << EOF > NOTS_LOWT.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=2


source ~/.bashrc

module purge

module load GCC/11.3.0  OpenMPI/4.1.4 PLUMED/2.8.1 CUDA/11.7.0
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi_nogpu -i $1/in.$3-out-lowT >> OUT 2>&1
EOF
}

Prepare_NOTS_HIGHT () {
		
cat << EOF > NOTS_HIGHT.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=2

#SBATCH -d afterany:$4

source ~/.bashrc

module purge

module load GCC/11.3.0  OpenMPI/4.1.4 PLUMED/2.8.1 CUDA/11.7.0 
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi_nogpu -i $1/in.$3-out-highT >> OUT 2>&1
EOF
}

Prepare_NOTS_FLOW () {
		
cat << EOF > NOTS_FLOW.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=2

#SBATCH -d afterany:$4
source ~/.bashrc

module purge

module load GCC/11.3.0  OpenMPI/4.1.4 PLUMED/2.8.1 CUDA/11.7.0 
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi_nogpu -i $1/in.$3-srd  >> OUT 2>&1
EOF
}



Prepare_KAROLINA_LOWT () {
		
cat << EOF > KAROLINA_LOWT.sh
#!/bin/bash
#SBATCH --job-name $2
#SBATCH --account DD-23-173
#SBATCH --partition qgpu
#SBATCH --time 24:00:00
#SBATCH --nodes=1
#SBATCH --gpus 4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20




module purge
module load  GROMACS/2021.4-fosscuda-2020b-PLUMED-2.7.3

#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 20 lmp_mpi  -sf gpu -pk gpu 4 neigh no -i ./in.$3-out-lowT >> OUT 2>&1
EOF
}


Prepare_KAROLINA_HIGHT () {
		
cat << EOF > KAROLINA_HIGHT.sh
#!/bin/bash
#SBATCH --job-name $2
#SBATCH --account DD-23-173
#SBATCH --partition qgpu
#SBATCH --time 24:00:00
#SBATCH --nodes=1
#SBATCH --gpus 4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH -d afterany:$4



module purge
module load  GROMACS/2021.4-fosscuda-2020b-PLUMED-2.7.3

#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 20 lmp_mpi  -sf gpu -pk gpu 4 neigh no -i ./in.$3-out-highT >> OUT 2>&1
EOF
}

Prepare_KAROLINA_FLOW () {
		
cat << EOF > KAROLINA_FLOW.sh
#!/bin/bash
#SBATCH --job-name $2
#SBATCH --account DD-23-173
#SBATCH --partition qgpu
#SBATCH --time 24:00:00
#SBATCH --gpus 4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH -d afterany:$4



module purge
module load  GROMACS/2021.4-fosscuda-2020b-PLUMED-2.7.3

#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 20 lmp_mpi  -sf gpu -pk gpu 4 neigh no -i ./in.$3-srd  >> OUT 2>&1
EOF
}



#ARGUMENTS     #1       #2        #3
#          PDB_FILE PARAM_FILE PKA_FILE

pdb=$(realpath $1)
param=$(realpath $2)
pka=$(realpath $3)
HPC=$4
MASK_LINKER=$(realpath $5)

proteina=$(echo $pdb|xargs -n 1 basename | awk -F ".pdb" '{print $1}')
geometry=$(grep Geometry $param | awk '{print $2}')
type=$(grep Simul_Type $param | awk '{print $2}')

NAME=$(echo $proteina-$geometry-$type)



mkdir -p $type-$geometry/$proteina
cd $type-$geometry/$proteina
start=$PWD

latest=$(echo Simul_`date +%d-%m_%H%M%S`)

for seed in {1..2}
do
for tempera in {2..10..4}
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
cp $MASK_LINKER ./mask_linker.dat
cp $pka ./
cp $param ./param_input.dat

##################### SET PARAMS #############################
sed -i "s+Seed=.*#+Seed= $seed #+g" param_input.dat
sed -i "s+Temp=.*#+Temp= $tempera #+g" param_input.dat
sed -i "s+Polymer_Fraction=.*#+Polymer_Fraction= $fraction #+g" param_input.dat
################### CREATE SIMULATION FILES ##################
GENERAL_PREP.sh $pdb param_input.dat $pka $ph
################### SUBMIT JOB #######################

ARG4_UPPER=$(echo "$HPC" | tr '[:lower:]' '[:upper:]')

if [[ "$ARG4_UPPER" == "HYPERION" ]]; then
    echo "HYPERION"
    Prepare_HYPERION_LOWT "$place" "$NAME" "$proteina"
    sbatch HYPERION_LOWT.sh > jobid_lowT.dat
    job=$(awk '{print $NF}' jobid_lowT.dat)
		Prepare_HYPERION_HIGHT "$place" "$NAME" "$proteina" "$job" 
    sbatch HYPERION_HIGHT.sh > jobid_highT.dat
    job=$(awk '{print $NF}' jobid_highT.dat)
    Prepare_HYPERION_FLOW "$place" "$NAME" "$proteina" "$job"
    sbatch HYPERION_FLOW.sh > jobid_flow.dat
fi
if [[ "$ARG4_UPPER" == "KAROLINA" ]]; then
    echo "KAROLINA"
    Prepare_KAROLINA_LOWT "$place" "$NAME" "$proteina"
    sbatch KAROLINA_LOWT.sh > jobid_lowT.dat
    job=$(awk '{print $NF}' jobid_lowT.dat)
		Prepare_KAROLINA_HIGHT "$place" "$NAME" "$proteina" "$job"
    sbatch KAROLINA_HIGHT.sh > jobid_highT.dat
    job=$(awk '{print $NF}' jobid_highT.dat)
    Prepare_KAROLINA_FLOW "$place" "$NAME" "$proteina" "$job"
    sbatch KAROLINA_FLOW.sh > jobid_flow.dat
fi
if [[ "$ARG4_UPPER" == "NOTS" ]]; then
    echo "NOTS"
    Prepare_NOTS_LOWT "$place" "$NAME" "$proteina"
    sbatch NOTS_LOWT.sh > jobid_lowT.dat
    job=$(awk '{print $NF}' jobid_lowT.dat)
		Prepare_NOTS_HIGHT "$place" "$NAME" "$proteina" "$job"
    sbatch NOTS_HIGHT.sh > jobid_highT.dat
    job=$(awk '{print $NF}' jobid_highT.dat)
    Prepare_NOTS_FLOW "$place" "$NAME" "$proteina" "$job"
    sbatch NOTS_FLOW.sh > jobid_flow.dat
fi

cd $start
done
done
done
done




