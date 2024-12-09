#!/bin/bash
#
#SBATCH --nodelist=node1,node4
#SBATCH --time=24:00:00
#SBATCH --partition day-long
# Outputs ----------------------------------
#SBATCH -o /home/%u/log/%x-%j.out
#SBATCH -e /home/%u/log/%x-%j.err
# ------------------------------------------

cd /home/data/eccolab/SPLaT_fMRI
# edit the R code in the Rscript string below depending on what targets you need to run
Rscript -e "targets::tar_make(starts_with('level2'),use_crew=FALSE)"
