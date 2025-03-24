#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 PDB_FILE PARAM_FILE SUBMISSION_SCRIPT"
    exit 1
fi

# Assign arguments to variables
PDB_FILE=$1
PARAM_FILE=$2
SUBMISSION_SCRIPT=$3

# Function to check if the submission script contains PBS options
is_pbs() {
    grep -q '#PBS' "$1"
}

# Function to check if the submission script contains SLURM options
is_slurm() {
    grep -q '#SBATCH' "$1"
}

# Function to submit the job using the appropriate command
submit_job() {
    if is_pbs "$1"; then
        qsub "$1"
    elif is_slurm "$1"; then
        sbatch "$1"
    else
        echo "Error: Unrecognized submission script format. Please use PBS (#PBS) or SLURM (#SBATCH) options."
        exit 1
    fi
}

tipo=$(echo $PDB_FILE | awk -F ".pdb" '{print $1}')
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
            cp $PDB_FILE $simul
            cp $PARAM_FILE $simul/param_input.dat

            NRESIDUE=$(awk '$1=="ATOM"&&$3=="CA"{print $6}' $PDB_FILE | tail -n 1)

            awk -v nres=$NRESIDUE 'NR>1&&$1<nres{m=m+1}END{print m}' mask.dat >  $simul/mask.dat
            awk -v nres=$NRESIDUE 'NR>1&&$1<nres{print $1}' mask.dat >> $simul/mask.dat

            awk -v nres=$NRESIDUE 'NR>1&&$1<nres{m=m+1}END{print m}' mask2.dat >  $simul/mask2.dat
            awk -v nres=$NRESIDUE 'NR>1&&$1<nres{print $1}' mask2.dat >> $simul/mask2.dat

            cp $SUBMISSION_SCRIPT $simul
            echo $seed
            cd $simul
            place=$(pwd)
            echo $place

            ##################### SET PARAMS #############################
            sed -i "s+Seed=.*#+Seed= $seed #+g" param_input.dat
            sed -i "s+Temp=.*#+Temp= $tempera #+g" param_input.dat
            sed -i "s+Polymer_Fraction=.*#+Polymer_Fraction= $fraction #+g" param_input.dat
            sed -i "s+FOLDER+$place+g" $SUBMISSION_SCRIPT
            sed -i "s+INFILE+in.${tipo}+g" $SUBMISSION_SCRIPT
            sed -i "s+NAME+${tipo}-P-${seed}-T-${tempera}-T-${fraction}+g" $SUBMISSION_SCRIPT
            ################### CREATE SIMULATION FILES ##################
            GRAFTFLOW-BRUSH-COMPLEX.sh $PDB_FILE param_input.dat
            ################### SUBMIT JOB #######################
            submit_job $SUBMISSION_SCRIPT > jobid.dat
            cd $start
        done
    done
done
