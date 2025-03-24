#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=NAME
#SBATCH --partition=regular
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=2


source ~/.bashrc

module purge
module load OpenMPI/3.1.4-GCC-8.3.0
module load GSL/2.6-GCC-8.3.0
module load Python/2.7.16-GCCcore-8.3.0
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd FOLDER
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi -i INFILE >> OUT 2>&1
