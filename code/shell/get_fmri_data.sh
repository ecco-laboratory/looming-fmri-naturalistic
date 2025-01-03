#!/bin/bash
read -p "Enter subj num(s) to pull (no leading 0s yes spaces): " SUBJ_NUMS

for SUBJ_NUM in ${SUBJ_NUMS[@]}; do
    printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM

    rsync -r -v mthieu@csic.som.emory.edu:/home/mthieu/kragel_data/FERN-KRAGELSPLAT/splat_${SUBJ_NUM_4D}* /archival/projects/SPLaT/data/fmri/dicom/
done
