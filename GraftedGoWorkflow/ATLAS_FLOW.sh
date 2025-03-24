#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=NAME
#SBATCH --partition=regular
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH --gres=gpu:p40:2
#SBATCH --mem=300gb
PREVJOB

source ~/.bashrc

module purge
module load OpenMPI/4.0.5-gcccuda-2020b    GROMACS/2021.3-fosscuda-2020b-PLUMED-2.7.2
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd FOLDER
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 20 lmp_mpi-p40-dev  -sf gpu -pk gpu 1 neigh no -i INFILE >> OUT 2>&1
