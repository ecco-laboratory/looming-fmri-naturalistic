#!/bin/bash
read -p "Enter your csic.som.emory.edu username: " CSIC_USERNAME

# This script is currently set up to FULLY SYNC THE CSIC AND ARCHIVAL VERSIONS OF THE FOLDER
# so if there is any data in the CSIC folder you don't want to download... DEAL WITH IT! 
# it may take a while to index the files at the beginning but then the --ignore-existing flag should prompt it to download only new subjects' data
# so it won't re-download the entire project's worth of dicoms every time
rsync -av --ignore-existing ${CSIC_USERNAME}@csic.som.emory.edu:/home/${CSIC_USERNAME}/kragel_data/FERN-KRAGELSPLAT/ /archival/projects/SPLaT/data/fmri/dicom/
