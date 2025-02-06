#!/bin/bash
#
#SBATCH --account=default
#SBATCH --exclude=gpu1,gpu2
#SBATCH --time=7-00:00:00
#SBATCH --mem=1M
#SBATCH --partition week-long  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o /home/%u/log/%x-%A.out
#SBATCH -e /home/%u/log/%x-%A.err
# ------------------------------------------
# take the job number to monitor as an argument
# the log file delimits the array number with a dash
# but the job id is with an underscore
# so for the log file, replace underscore with dash
echo "Monitoring job: $1"
JOB_ID_LOG=$(echo $1 | tr _ -)
ARRAY_ID=$(echo $1 | cut -d '_' -f 2)
LOG_FILE_PATH="/home/mthieu/log/sbatch_fmriprep.sh-${JOB_ID_LOG}.err"
while true
do
    # find the error message in the log file
    # if the error message is present, cancel the job and restart both fmriprep and the monitoring process
    # error message is present if grep output has more than 0 lines
    GREP_LENGTH=$(grep "concurrent.futures.process.BrokenProcessPool" $LOG_FILE_PATH | wc -l)
    if [ $GREP_LENGTH -gt 0 ]; then 
        scancel $1
        echo "error found, cancelled job"
        NEW_SBATCH_OUT=$(sbatch --array=$ARRAY_ID /home/data/eccolab/SPLaT_fMRI/code/shell/sbatch_fmriprep.sh)
        NEW_JOB_ID=${NEW_SBATCH_OUT: -5}_${ARRAY_ID}
        sbatch /home/data/eccolab/SPLaT_fMRI/code/shell/monitor_fmriprep.sh $NEW_JOB_ID
        exit
    fi
    # ideally: also monitor for whether the main fmriprep job has finished
    # if it has finished, end this script as well (and don't restart it)
    if [ $(squeue | grep $1 | wc -l) -eq 0 ]; then 
        echo "it finished"
        exit
    fi
    # repeat every 5 seconds until script or job has ended
    echo "waiting..."
    sleep 5
done
