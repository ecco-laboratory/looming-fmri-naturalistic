#!/bin/bash
#
#SBATCH --account=default
#SBATCH --exclude=node3,gpu2
#SBATCH --time=7-00:00:00
#SBATCH --mem=48G
#SBATCH --partition week-long  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o /home/%u/log/%x-%A-%a.out
#SBATCH -e /home/%u/log/%x-%A-%a.err
#SBATCH --mail-user=mthieu@emory.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

STUDY_DIR="/home/data/eccolab/SPLaT_fMRI/"
BIDS_DIR="${STUDY_DIR}ignore/data/fmri/"
DERIVS_DIR="derivatives/fmriprep-23.1.4"

# Prepare some writeable bind-mount points.
TEMPLATEFLOW_HOST_HOME=$HOME/.cache/templateflow
FMRIPREP_HOST_CACHE=$HOME/.cache/fmriprep
mkdir -p ${TEMPLATEFLOW_HOST_HOME}
mkdir -p ${FMRIPREP_HOST_CACHE}

# Prepare derivatives folder
mkdir -p ${BIDS_DIR}/${DERIVS_DIR}


# Designate a templateflow bind-mount point
export SINGULARITYENV_TEMPLATEFLOW_HOME="/templateflow"
SINGULARITY_CMD="singularity run --cleanenv --bind ${BIDS_DIR}:/data --bind /home/data/shared/SingularityImages/:/fslicensepath --bind ${TEMPLATEFLOW_HOST_HOME}:${SINGULARITYENV_TEMPLATEFLOW_HOME}  --bind ${STUDY_DIR}/ignore/tmp:/work  /home/data/shared/SingularityImages/fmriprep-23.1.4.simg"

echo Current slurm array task ID: $SLURM_ARRAY_TASK_ID
# Parse the participants.tsv file and extract one subject ID from the line corresponding to this SLURM task.
subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )

echo Running for subject: $subject
# Remove IsRunning files from FreeSurfer
# find ${BIDS_DIR}/derivatives/freesurfer-6.0.1/sub-$subject/ -name "*IsRunning*" -type f -delete

# Compose the command line
# flags you may or may not want to use:
# --fs-no-reconall (default would be to run it)
# flag notes
# --nprocs is the total number of threads
# --omp-nthreads is the total number of threads per process
# I think we don't use the --mem flag because that is capped by the memory requested in the sbatch header?
cmd="${SINGULARITY_CMD} /data /data/${DERIVS_DIR} participant --participant-label $subject --fs-license-file /fslicensepath/license.txt --fs-no-reconall -w /work/ -vv --nprocs 32 --omp-nthreads 8  --output-spaces MNI152NLin2009cAsym:res-2 anat fsaverage5 --verbose --low-mem --me-output-echos"


# Setup done, run the command
echo Running task ${SLURM_ARRAY_TASK_ID}
echo Commandline: $cmd
eval $cmd
exitcode=$?

# Output results to a table
echo "sub-$subject   ${SLURM_ARRAY_TASK_ID}    $exitcode" \
      >> ${SLURM_JOB_NAME}.${SLURM_ARRAY_JOB_ID}.tsv
echo Finished tasks ${SLURM_ARRAY_TASK_ID} with exit code $exitcode
exit $exitcode
