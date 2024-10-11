#!/bin/bash
# It only runs FROM Monica's local TO server rn

read -p "Enter subj num(s) to heuristic with no leading 0s (but with spaces): " SUBJ_NUMS

for SUBJ_NUM in ${SUBJ_NUMS[@]}; do
    printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM
# can't scp into symlinked folders, must go to the original
    scp ~/OneDrive\ -\ Emory\ University/PsychoPy/looming_nback_psy/data/sub-${SUBJ_NUM_4D}/*.csv ecco:/archival/projects/SPLaT/data/beh/sub-${SUBJ_NUM_4D}/raw/
    scp ~/OneDrive\ -\ Emory\ University/PsychoPy/looming_stimrating_psy/data/sub-${SUBJ_NUM_4D}/*.csv ecco:/archival/projects/SPLaT/data/beh/sub-${SUBJ_NUM_4D}/raw/
done
