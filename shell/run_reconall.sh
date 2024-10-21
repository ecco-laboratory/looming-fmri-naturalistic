#!/bin/zsh
# It only runs on Monica's local (where freesurfer is installed)

read -p "Enter subj num to freesurf with no leading 0s: " SUBJ_NUM
printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM

scp ecco:/archival/projects/SPLaT/data/fmri/nifti/sub-${SUBJ_NUM_4D}/anat/sub-${SUBJ_NUM_4D}_T1w.nii.gz /Applications/freesurfer/7.4.1/subjects/
# this one takes hours to run so you may never be able to run this whole script as .sh
# you may just have to copy the commands into console as you need them
recon-all -s sub-${SUBJ_NUM_4D} -i $SUBJECTS_DIR/sub-${SUBJ_NUM_4D}_T1w.nii.gz -all
mris_convert $SUBJECTS_DIR/sub-${SUBJ_NUM_4D}/surf/rh.pial $SUBJECTS_DIR/sub-${SUBJ_NUM_4D}/surf/rh.stl
mris_convert $SUBJECTS_DIR/sub-${SUBJ_NUM_4D}/surf/lh.pial $SUBJECTS_DIR/sub-${SUBJ_NUM_4D}/surf/lh.stl
