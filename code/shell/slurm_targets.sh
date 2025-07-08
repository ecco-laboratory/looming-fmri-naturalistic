#!/bin/bash
#
#SBATCH --account=default
#SBATCH --exclude=gpu1,gpu2
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --partition day-long
# Outputs ----------------------------------
#SBATCH -o /home/%u/log/%x-%j.out
#SBATCH -e /home/%u/log/%x-%j.err
# ------------------------------------------

# takes the TAR_PROJECT as a positional arg now!
cd /home/data/eccolab/SPLaT_fMRI
# edit the R code in the Rscript string below depending on what targets you need to run
# random notes:
# the targets that call canlabtools preproc (bold.masked and wb.model.connectivity by-subject targets) require like 32 GB memory
Rscript -e "Sys.setenv(TAR_PROJECT='${1}'); library(targets); tar_make(contains('${2}'), use_crew=FALSE)"
