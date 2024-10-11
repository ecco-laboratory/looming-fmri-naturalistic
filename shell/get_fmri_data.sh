#!/bin/bash
read -p "Enter subj number to pull (with leading 0s): " SUBJ_NUM
# TODO: Wildcard it so it will pull the subject but we don't need to specify the whole long folder name
rsync -r mthieu@csic.som.emory.edu:/home/mthieu/kragel_data/FERN-KRAGELSPLAT/splat_${SUBJ_NUM}* /archival/projects/SPLaT/data/fmri/dicom/
