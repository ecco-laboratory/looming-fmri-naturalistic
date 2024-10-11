#!/bin/bash
# TODO: decide whether you are going to use a standard heuristic file or a new one for each sub like CVA

read -p "Enter subj num(s) to heuristic with no leading 0s (but with spaces): " SUBJ_NUMS

for SUBJ_NUM in ${SUBJ_NUMS[@]}; do
    printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM
    singularity run -B /archival/projects/SPLaT/data/fmri:/base \
    /home/data/shared/SingularityImages/heudiconv_1.0.0.sif \
    -d /base/dicom/splat_{subject}*/*/*.dcm \
    -o /base/nifti \
    -f /base/heuristic_default.py \
    -s $SUBJ_NUM_4D \
    -c dcm2niix \
    -b \
    --overwrite
done
