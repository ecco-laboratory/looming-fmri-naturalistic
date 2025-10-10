#!/bin/bash
#
#SBATCH --account=default
#SBATCH --time=0-24:00:00
#SBATCH --mem=24G
#SBATCH --partition day-long  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o /home/%u/log/%x-%A-%a.out
#SBATCH -e /home/%u/log/%x-%A-%a.err
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

singularity run -B /archival/projects/SPLaT/data/mesh:/base \
/home/data/shared/SingularityImages/freesurfer_8.0.0_20250210.simg \
mris_convert /base/sub-${SUBJ_NUM_4D}/surf/rh.pial /base/sub-${SUBJ_NUM_4D}/surf/rh.stl

singularity run -B /archival/projects/SPLaT/data/mesh:/base \
/home/data/shared/SingularityImages/freesurfer_8.0.0_20250210.simg \
mris_convert /base/sub-${SUBJ_NUM_4D}/surf/lh.pial /base/sub-${SUBJ_NUM_4D}/surf/lh.stl

/home/data/eccolab/CondaEnvs/pymeshlab/bin/python /home/data/eccolab/SPLaT_fMRI/code/python/combine_downsample_recon_mesh.py --subj_num $SLURM_ARRAY_TASK_ID
