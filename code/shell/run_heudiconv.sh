#!/bin/bash
# This study uses a single heuristic file that should be consistent across all subjects
# because they all get scanned with the exact same sequence
# This also only runs in an interactive bash instance, not in slurm
# Monica tried setting it up to run through slurm and it appeared to convert the data okay
# but then didn't update the participants.tsv file with the finished subjects??

read -p "Enter subj num to heuristic with no leading 0s (but with spaces): " SUBJ_NUM

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

