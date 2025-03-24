#!/bin/bash
#PBS -q qgpu_eurohpc
#PBS -N NAME
#PBS -l select=1:mpiprocs=128:ngpus=8 
#PBS -A DD-23-101
#PBS -l walltime=24:00:00

 

SCRDIR=FOLDER
cd $SCRDIR || exit

 
module purge
module load  GROMACS/2021.4-fosscuda-2020b-PLUMED-2.7.3
module list
nvidia-smi  >> OUT 2>&1 
mpirun -np 20 lmp_mpi -sf gpu -pk gpu 1 neigh no -i INFILE >> OUT 2>&1
exit
