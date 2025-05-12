#!/bin/bash
# Now only downloads server-pymeshlab'd FINISHED mesh to local for emailing to participant

read -p "Enter subj num to download with no leading 0s: " SUBJ_NUM
printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM

scp ecco:/archival/projects/SPLaT/data/mesh/sub-${SUBJ_NUM_4D}/surf/bl.stl /Users/mthieu/Downloads/sub-${SUBJ_NUM_4D}.stl
