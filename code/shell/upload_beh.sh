#!/bin/bash
# It only runs FROM Elias' local TO server rn

read -p "Enter subj num(s) to upload with no leading 0s (but with spaces): " SUBJ_NUMS

for SUBJ_NUM in ${SUBJ_NUMS[@]}; do
    printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM
# can't scp into symlinked folders, must go to the original
    scp ~/OneDrive\ -\ Emory\ University/SPLaT/data/controlled/sub-${SUBJ_NUM_4D}/*.csv ecco:/archival/projects/SPLaT/data/beh/sub-${SUBJ_NUM_4D}/raw/
    scp ~/OneDrive\ -\ Emory\ University/SPLaT/data/naturalistic/sub-${SUBJ_NUM_4D}/*.csv ecco:/archival/projects/SPLaT/data/beh/sub-${SUBJ_NUM_4D}/raw/
done
