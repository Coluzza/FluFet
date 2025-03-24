#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=NAME
#SBATCH --partition=gpu
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH --gres=gpu:4

source ~/.bashrc

module purge
module load GROMACS/2020.3-foss-2021a-PLUMED-2.7.2 
module load OpenMPI/4.0.5-gcccuda-2020b
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd FOLDER
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 20 lmp_mpi -sf gpu -i INFILE >> OUT 2>&1
