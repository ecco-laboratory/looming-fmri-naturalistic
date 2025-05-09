#!/bin/bash
#
#SBATCH --account=default
#SBATCH --nodelist=gpu1
#SBATCH --time=0-24:00:00
#SBATCH --mem=24G
#SBATCH --gres=gpu:1
#SBATCH --partition day-long  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o /home/%u/log/%x-%A.out
#SBATCH -e /home/%u/log/%x-%A.err
#SBATCH --mail-user=mthieu@emory.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

# This sbatch script takes the subj nums to freesurf, with no leading 0s, 
# as the sbatch --array argument 
printf -v SUBJ_NUM_4D '%04d' $SLURM_ARRAY_TASK_ID

singularity run -B /archival/projects/SPLaT/data:/base \
--env SUBJECTS_DIR=/base/mesh,FS_ALLOW_DEEP=1 \
/home/data/shared/SingularityImages/freesurfer_8.0.0_20250210.simg \
recon-all \
-s sub-${SUBJ_NUM_4D} \
-i /base/fmri/nifti/sub-${SUBJ_NUM_4D}/anat/sub-${SUBJ_NUM_4D}_T1w.nii.gz \
-all
