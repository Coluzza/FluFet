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
#SBATCH --gres=gpu:rtx3090:1

source ~/.bashrc

module purge

module load LAMMPS/23Jun2022-foss-2021b-kokkos-CUDA-11.4.1
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.$3-out-lowT >> OUT 2>&1
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
#SBATCH --gres=gpu:rtx3090:1
#SBATCH -d afterany:$4

source ~/.bashrc

module purge

module load LAMMPS/23Jun2022-foss-2021b-kokkos-CUDA-11.4.1
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 2 lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.$3-out-highT >> OUT 2>&1
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
#SBATCH --gres=gpu:rtx3090:1
#SBATCH -d afterany:$4
source ~/.bashrc

module purge

module load LAMMPS/23Jun2022-foss-2021b-kokkos-CUDA-11.4.1
#CHANGE FOLDER TO THE FULL PATH OF YOUR SIMULATION FOLDER
cd $1
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS
mpirun -np 10 lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.$3-srd  >> OUT 2>&1
EOF
}

Prepare_NOTS_MIN () {
cat << EOF > NOTS_MIN.sh
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
srun lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.$3-out-min >> OUT 2>&1
EOF
}




Prepare_NOTS_LOWT () {
cat << EOF > NOTS_LOWT.sh
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
srun lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.$3-out-lowT >> OUT 2>&1
EOF
}

Prepare_NOTS_HIGHT () {
		
cat << EOF > NOTS_HIGHT.sh
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
srun lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.$3-out-highT >> OUT 2>&1
EOF
}

Prepare_NOTS_FLOW () {
		
cat << EOF > NOTS_FLOW.sh
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
srun lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.$3-srd  >> OUT 2>&1
EOF
}

Prepare_NOTS_FLOW_RESTART () {
		
cat << EOF > NOTS_FLOW_RESTART.sh
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
srun lmp_mpi  -sf gpu -pk gpu 1 neigh no -i $1/in.$3-srd-restart  >> OUT 2>&1
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
#=========================================== END OLD SCRIPTS =====================================



Prepare_NOTS_FIRST () {
cat << EOF > NOTS_FIRST.sh
#!/bin/bash

#SBATCH --time=24:00:00  # Wall time
#SBATCH --job-name=$2
#SBATCH --account=commons
#SBATCH --partition=commons
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
#SBATCH --account=commons
#SBATCH --partition=commons
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


Prepare_COGITATORE_FIRST () {
cat << EOF > COGITATORE_FIRST.sh
#!/bin/bash

#SBATCH --job-name=$2
#SBATCH --error=Standard_%j.err  # stderr
#SBATCH --output=Standard_%j.out  # stdout
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=4                      # Number of tasks (single task)
source ~/.bashrc


mkdir -p /data_storage/scratch/
# Define paths
MASTER_SIM_DIR="$1"  # Simulation folder on the master
#SLAVE_SIM_DIR=\$(mktemp -d /data_storage/scratch/"$2"_XXXXX)
SLAVE_SIM_DIR="/data_storage/scratch/\$SLURM_JOB_ID"

echo MASTER_SIM_DIR \$MASTER_SIM_DIR
echo SLAVE_SIM_DIR \$SLAVE_SIM_DIR

MASTER_NODE="Cogitatore1"  # Master node hostname or IP

# Copy the simulation folder from the master to the compute node
echo "Copying simulation folder to compute node..."

mkdir -p \${SLAVE_SIM_DIR}/
echo \$MASTER_SIM_DIR > \${SLAVE_SIM_DIR}/master_dir.txt
# Check if the directory was created successfully and is not root or home
if [[ ! -d \${SLAVE_SIM_DIR} || \${SLAVE_SIM_DIR} == "/" || \${SLAVE_SIM_DIR} == "\$HOME" ]]; then
	echo "Error: Failed to create slave directory or directory is set to root/home."
	echo SLAVE_SIM_DIR \$SLAVE_SIM_DIR
	exit 1
fi
rsync -auvC cogitatore1:\${MASTER_SIM_DIR}/ \${SLAVE_SIM_DIR}/


# Regularly sync the results back to the master node
echo "Syncing results back to master node..."

while true; do
	rsync -auvC \${SLAVE_SIM_DIR}/ cogitatore1:\${MASTER_SIM_DIR}/
	sleep 20m  # Sync every 20 minutes (adjust as needed)
done &


# Run the simulation (replace this with your actual simulation command)
echo "Starting simulation on \$SLURMD_NODENAME..."
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS

cd \${SLAVE_SIM_DIR}/

mpirun -np 4 lmp_mpi  -sf gpu -pk gpu 1 neigh no -i \${SLAVE_SIM_DIR}/$3 >> OUT 2>&1

# Final sync after the simulation completes
echo "Final sync to master node..."

rsync -auvC  --remove-source-files \${SLAVE_SIM_DIR}/ cogitatore1:\${MASTER_SIM_DIR}/

echo "Simulation completed on \$SLURMD_NODENAME."

EOF
}

Prepare_COGITATORE_NEXT () {
cat << EOF > COGITATORE_NEXT.sh
#!/bin/bash

#SBATCH --job-name=$2
#SBATCH --error=Standard_%j.err  # stderr
#SBATCH --output=Standard_%j.out  # stdout
#SBATCH --nodes=1
#SBATCH --ntasks=4                      # Number of tasks (single task)
#SBATCH -d afterany:$4

source ~/.bashrc

# Define paths
MASTER_SIM_DIR="$1"  # Simulation folder on the master
#SLAVE_SIM_DIR=\$(mktemp -d /data_storage/scratch/"$2"_XXXXX)
SLAVE_SIM_DIR="/data_storage/scratch/\$SLURM_JOB_ID"

echo MASTER_SIM_DIR \$MASTER_SIM_DIR
echo SLAVE_SIM_DIR \$SLAVE_SIM_DIR


MASTER_NODE="Cogitatore1"                        # Master node hostname or IP

# Copy the simulation folder from the master to the compute node
echo "Copying simulation folder to compute node..."

mkdir -p \${SLAVE_SIM_DIR}/
echo \$MASTER_SIM_DIR > \${SLAVE_SIM_DIR}/master_dir.txt
# Check if the directory was created successfully and is not root or home
if [[ ! -d \${SLAVE_SIM_DIR} || \${SLAVE_SIM_DIR} == "/" || \${SLAVE_SIM_DIR} == "\$HOME" ]]; then
	echo "Error: Failed to create slave directory or directory is set to root/home."
	echo SLAVE_SIM_DIR \$SLAVE_SIM_DIR
	exit 1
fi
rsync -auvC cogitatore1:\${MASTER_SIM_DIR}/ \${SLAVE_SIM_DIR}/


# Regularly sync the results back to the master node
echo "Syncing results back to master node..."

while true; do
	rsync -auvC \${SLAVE_SIM_DIR}/ cogitatore1:\${MASTER_SIM_DIR}/
	sleep 20m  # Sync every 20 minutes (adjust as needed)
done &


# Run the simulation (replace this with your actual simulation command)
echo "Starting simulation on \$SLURMD_NODENAME..."
#CHANGE INFILE TO THE NAME OF THE LAMMPS in. FILE WITH THE SIMUALTION PARAMETERS

cd \${SLAVE_SIM_DIR}/

mpirun -np 4 lmp_mpi  -sf gpu -pk gpu 1 neigh no -i \${SLAVE_SIM_DIR}/$3 >> OUT 2>&1

# Final sync after the simulation completes
echo "Final sync to master node..."

rsync -auvC --remove-source-files \${SLAVE_SIM_DIR}/ cogitatore1:\${MASTER_SIM_DIR}/
echo "Simulation completed on \$SLURMD_NODENAME."
EOF
}

Prepare_EPILOG () {
cat << EOF > EPILOG.sh

#!/bin/bash

LOGFILE="/var/log/slurm/epilog.log"

echo "Epilog started for JobID=\$SLURM_JOB_ID on Node=\$SLURMD_NODENAME" >> "\$LOGFILE"

# Retrieve the simulation directories
MASTER_SIM_DIR="\$SLURM_SUBMIT_DIR"  # Job submission directory
SLAVE_SIM_DIR="/data_storage/scratch/\$SLURM_JOB_ID"


# **Process Cleanup Section**
# Get the Slurm step daemon PID
pids=\$(pgrep -f "\$SLURM_JOB_ID")

    

if [[ -n "\$pids" ]]; then
    echo "Killing processes for JobID=\$SLURM_JOB_ID: \$pids" >> "\$LOGFILE"
    kill -9 \$pids 2>/dev/null
else
    echo "No running processes found for JobID=\$SLURM_JOB_ID" >> "\$LOGFILE"
fi




echo "MASTER_SIM_DIR: \$MASTER_SIM_DIR" >> "\$LOGFILE"
echo "SLAVE_SIM_DIR: \$SLAVE_SIM_DIR" >> "\$LOGFILE"

echo "Syncing simulation results from \$SLAVE_SIM_DIR to master at \$MASTER_SIM_DIR" >> "\$LOGFILE"
# Ensure directories exist before syncing

if [[ ! -d \${SLAVE_SIM_DIR} || \${SLAVE_SIM_DIR} == "/" || \${SLAVE_SIM_DIR} == "\$HOME" ]]; then
	echo "Error: Slave directory does not exists or directory is set to root/home."
	echo SLAVE_SIM_DIR \$SLAVE_SIM_DIR
	exit 1
fi


# Rsync simulation results and remove source files after successful transfer
rsync -auvC --remove-source-files "\${SLAVE_SIM_DIR}/" "cogitatore1:\${MASTER_SIM_DIR}/"

# Remove empty directories after rsync
rm -rf "\$SLAVE_SIM_DIR"

echo "Final sync to master node completed." >> "\$LOGFILE"

echo "Epilog completed for JobID=\$SLURM_JOB_ID" >> "\$LOGFILE"


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
geometry=$(grep Geometry $param  | awk '{gsub(/Geometry\s*=/,"",$0);print $1}')
type=$(grep Simul_Type $param | awk '{gsub(/Simul_Type\s*=/,"",$0);print $1}')

NAME=$(echo $proteina-$geometry-$type)



mkdir -p $type-$geometry/$proteina
cd $type-$geometry/$proteina
start=$PWD

latest=$(echo Simul_`date +%d-%m_%H%M%S`)

mkdir $latest
echo $@ > $latest/logfile.dat

# Read parameters from $param
seeds=$(grep "Seed" $param | awk '{gsub(/Seed\s*=/,"",$0); gsub(/#.*/,"",$0); print $0}')
fractions=$(grep "Polymer_Fraction" $param | awk '{gsub(/Polymer_Fraction\s*=/,"",$0); gsub(/#.*/,"",$0); print $0}')
lengths=$(grep "Polymer_Length" $param | awk '{gsub(/Polymer_Length\s*=/,"",$0); gsub(/#.*/,"",$0); print $0}')
temps=$(grep "Final_Flow_Temp" $param | awk '{gsub(/Final_Flow_Temp\s*=/,"",$0); gsub(/#.*/,"",$0); print $0}')
escales=$(grep "E_Scale" $param | awk '{gsub(/E_Scale\s*=/,"",$0); gsub(/#.*/,"",$0); print $0}')
nprots=$(grep "N_Proteins" $param | awk '{gsub(/N_Proteins\s*=/,"",$0); gsub(/#.*/,"",$0); print $0}')
phs=$(grep "pH" $param | awk '{gsub(/pH\s*=/,"",$0); gsub(/#.*/,"",$0); print $0}')
pressures=$(grep "Pressure" $param | awk '{gsub(/Pressure\s*=/,"",$0); gsub(/#.*/,"",$0); print $0}')

# Convert space-separated values to arrays
IFS=' ' read -r -a seed_array <<< "$seeds"
IFS=' ' read -r -a fraction_array <<< "$fractions"
IFS=' ' read -r -a length_array <<< "$lengths"
IFS=' ' read -r -a temp_array <<< "$temps"
IFS=' ' read -r -a nprot_array <<< "$nprots"
IFS=' ' read -r -a ph_array <<< "$phs"
IFS=' ' read -r -a pressure_array <<< "$pressures"
IFS=' ' read -r -a escale_array <<< "$escales"

# Loop through each combination of parameters
for seed in "${seed_array[@]}"
do
for fraction in "${fraction_array[@]}"
do
for length in "${length_array[@]}"
do
for temp in "${temp_array[@]}"
do
for nprot in "${nprot_array[@]}"
do
for ph in "${ph_array[@]}"
do
for pressure in "${pressure_array[@]}"
do
for escale in "${escale_array[@]}"
do

simul=$latest"/P-"$seed"-F-"$fraction"-T-"$temp"-N-"$nprot"-pH-"$ph"-P-"$pressure-E-"$escale"-L-"$length"
mkdir -p $simul

cd $simul
place=$(pwd)
echo $place

cp $pdb ./
cp $MASK_LINKER ./mask_linker.dat
cp $pka ./
grep -v "pH" $param > ./param_input.dat

##################### SET PARAMS #############################
sed -i "s+Seed=.*#+Seed= $seed #+g" param_input.dat
sed -i "s+Final_Flow_Temp=.*#+Final_Flow_Temp= $temp #+g" param_input.dat
sed -i "s+N_Proteins=.*#+N_Proteins= $nprot #+g" param_input.dat
sed -i "s+Polymer_Fraction=.*#+Polymer_Fraction= $fraction #+g" param_input.dat
sed -i "s+Polymer_Length=.*#+Polymer_Length= $length #+g" param_input.dat
sed -i "s+Pressure=.*#+Pressure= $pressure #+g" param_input.dat
sed -i "s+E_Scale=.*#+E_Scale= $escale #+g" param_input.dat
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
		# Initialize counter
		i=1
		echo "NOTS"
		
		
		# Loop over files in numerical order from 1 to 10
for num in {1..17}; do
  # Construct the filename
  file="in.$proteina-$num"
  		# Check if the filename actually found files
  		if [[ -e $file ]]; then
    		echo "Processing $file"
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
    		# Increment counter
    		((i++))
  		fi
		done
fi

if [[ "$ARG4_UPPER" == "COGITATORE" ]]; then
	# Initialize counter
	i=1
	echo "COGITATORE"
	
	# Loop over files in numerical order from 1 to 10
	for num in {1..18}; do
 		# Construct the filename
  		file="in.$proteina-$num"
		# Check if the filename actually found files
		if [[ -e $file ]]; then
			echo "Processing $file"
				((k=i-1))
				if [[ $k -eq 0 ]]; then
					Prepare_COGITATORE_FIRST "$place" "$NAME" "$file" 
					mv COGITATORE_FIRST.sh COGITATORE_$i.sh
				
				else
					job=$(awk 'END{print $NF}'  jobid_$k.dat)
					Prepare_COGITATORE_NEXT "$place" "$NAME" "$file" "$job"
					mv COGITATORE_NEXT.sh COGITATORE_$i.sh
				fi
				sbatch COGITATORE_$i.sh > jobid_$i.dat
			# Increment counter
			((i++))
		fi
	done
	
fi

cd $start
done
done
done
done
done
done
done
done

exit
if [[ "$ARG4_UPPER" == "NOTS" ]]; then
    echo "NOTS"
		Prepare_NOTS_MIN "$place" "$NAME" "$proteina"
    sbatch NOTS_MIN.sh > jobid_min.dat
    job=$(awk '{print $NF}' jobid_min.dat)
		
    Prepare_NOTS_LOWT "$place" "$NAME" "$proteina" "$job"
    sbatch NOTS_LOWT.sh > jobid_lowT.dat
    job=$(awk '{print $NF}' jobid_lowT.dat)
		
		Prepare_NOTS_HIGHT "$place" "$NAME" "$proteina" "$job"
    sbatch NOTS_HIGHT.sh > jobid_highT.dat
    job=$(awk '{print $NF}' jobid_highT.dat)
		
    Prepare_NOTS_FLOW "$place" "$NAME" "$proteina" "$job"
    sbatch NOTS_FLOW.sh > jobid_flow-1.dat
		for j in {2..3}
		do
		((k=j-1))
		job=$(awk 'END{print $NF}'  jobid_flow-$k.dat)
		Prepare_NOTS_FLOW_RESTART "$place" "$NAME" "$proteina" "$job"
		mv NOTS_FLOW_RESTART.sh NOTS_FLOW_RESTART-$j.sh
		sbatch NOTS_FLOW_RESTART-$j.sh > jobid_flow-$j.dat
		done
		
		
fi

