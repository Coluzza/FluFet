# Set threshold from command line argument or default to 1
threshold=${1:-1}

get_last_file() {
	/bin/ls -d1 in.*  | grep -oP "in\..*-\d+" | sort -t- -k2,2n | tail -n 1
}
get_last_job() {
	/bin/ls -d1 jobid_*.dat | awk '{gsub("jobid_","",$0);gsub(".dat","",$0); print $0}' | sort -k 1n,1 | tail -n 1
}
Prepare_NOTS_FIRST () {
cat << EOF > NOTS_FIRST.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1


source ~/.bashrc

module purge

module load GCC/11.3.0  OpenMPI/4.1.4 PLUMED/2.8.1 CUDA/11.7.0
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
srun lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/$3 >> OUT 2>&1
EOF
}
Prepare_NOTS_NEXT () {
cat << EOF > NOTS_NEXT.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --error=Standard_%N.err  # stderr
#SBATCH --output=Standard_%N.out  # stdout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH -d afterany:$4

source ~/.bashrc

module purge

module load GCC/11.3.0  OpenMPI/4.1.4 PLUMED/2.8.1 CUDA/11.7.0
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
srun lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/$3 >> OUT 2>&1
EOF
}



echo "NOTS"
proteina=$(grep Prefix param.dat | awk '{gsub(/Prefix\s*=/,"",$0);print $1}')
place=$(pwd)
NAME=$(echo $proteina-$geometry-$type)

 if [[ $threshold -le 10 ]]; then

# Initialize counter
i=1

# Loop over files starting with "in-"
for num in {1..10}; do
  # Construct the filename
  file="in.$proteina-$num"  
# Check if the filename actually found files
  if [[ -e $file ]]; then
    echo "Processing $file $i"
    
    # Process only if i is greater than or equal to the threshold
    if [[ $i -ge $threshold ]]; then
      ((k=i-1))
      if [[ $k -eq 0 ]]; then
        Prepare_NOTS_FIRST "$place" "$NAME" "$file" 
        mv NOTS_FIRST.sh NOTS_$i.sh
      else
        job=$(awk 'END{print $NF}'  jobid_$k.dat)
        Prepare_NOTS_NEXT "$place" "$NAME" "$file" "$job"
        mv NOTS_NEXT.sh NOTS_$i.sh
      fi
      sbatch NOTS_$i.sh > jobid_$i.dat
    fi
    
    # Increment counter
    ((i++))
  fi
done
else
	k=$(get_last_job)

	((i=k+1))
	job=$(awk 'END{print $NF}'  jobid_$k.dat)
	Prepare_NOTS_NEXT "$place" "$NAME" "in.$proteina-10" "$job"
        mv NOTS_NEXT.sh NOTS_$i.sh
	sbatch "NOTS_$i.sh" > jobid_$i.dat
fi
