#!/bin/bash
LOGFILE="/var/log/slurm/epilog.log"
MASTER_SIM_DIR="SIMULATIONFOLDER"  # Simulation folder on the master
SLAVE_SIM_DIR=$(mktemp -d /data_storage/scratch/"$2"_XXXXX)

echo MASTER_SIM_DIR $MASTER_SIM_DIR
echo SLAVE_SIM_DIR $SLAVE_SIM_DIR

echo "Epilog started for JobID=$SLURM_JOB_ID on Node=$SLURMD_NODENAME" >> "$LOGFILE"

# Final sync after the simulation completes
echo "Final sync to master node..."

rsync -auvC --remove-source-files ${SLAVE_SIM_DIR}/ cogitatore1:${MASTER_SIM_DIR}/

# Check if the job ID exists
if [[ -n "$SLURM_JOB_ID" ]]; then
    # Get the list of PIDs related to this job
    #pids=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/JobId/{print $2}' | tr -d ' ')
    pids=$(pgrep -f "$SLURM_JOB_ID")

    # If PIDs are found, kill them
    if [[ -n "$pids" ]]; then
        echo "Killing processes for JobID=$SLURM_JOB_ID: $pids" >> "$LOGFILE"
        kill -9 $pids 2>/dev/null
    else
        echo "No running processes found for JobID=$SLURM_JOB_ID" >> "$LOGFILE"
    fi
else
    echo "No Job ID detected. Skipping process cleanup." >> "$LOGFILE"
fi

echo "Epilog completed for JobID=$SLURM_JOB_ID" >> "$LOGFILE"