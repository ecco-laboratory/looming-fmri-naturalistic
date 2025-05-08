#!/bin/bash
# this script takes the array ID that's going into sbatch_fmriprep as a positional argument
read -p "Enter subj num to fmriprep with no leading 0s: " SUBJ_NUM

SCRIPT_DIR=/home/data/eccolab/SPLaT_fMRI/code/shell

SBATCH_OUT_FMRIPREP=$(sbatch --array=$SUBJ_NUM $SCRIPT_DIR/sbatch_fmriprep.sh)
JOB_ID=${SBATCH_OUT_FMRIPREP: -5}_${SUBJ_NUM}
echo "${SBATCH_OUT_FMRIPREP}: main fmriprep job"
SBATCH_OUT_MONITOR=$(sbatch $SCRIPT_DIR/monitor_fmriprep.sh $JOB_ID)
echo "${SBATCH_OUT_MONITOR}: monitoring job"
