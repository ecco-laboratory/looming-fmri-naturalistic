#!/bin/bash
# Now only downloads server-pymeshlab'd FINISHED mesh to local for emailing to participant

read -p "Enter subj num to download with no leading 0s: " SUBJ_NUM
printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM

read -p "Enter path to download folder (abs or rel from curr local dir): " LOCAL_DL_FOLDER
scp ecco:/archival/projects/SPLaT/data/mesh/sub-${SUBJ_NUM_4D}/surf/bl.stl ${LOCAL_DL_FOLDER}/sub-${SUBJ_NUM_4D}.stl
