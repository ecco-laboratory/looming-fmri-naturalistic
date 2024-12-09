#!/bin/bash
read -p "Enter subj number to pull (no leading 0s): " SUBJ_NUM
printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM

rsync -r -v mthieu@csic.som.emory.edu:/home/mthieu/kragel_data/FERN-KRAGELSPLAT/splat_${SUBJ_NUM_4D}* /archival/projects/SPLaT/data/fmri/dicom/
