#!/bin/bash

# EDIT THIS TO MATCH YOUR LOCAL INSTALL PATH FOR OUR REPO!
REPO_PATH='/path/to/your/repo/install'
cd $REPO_PATH

# Make subdirectory structure matching my original setup so the relative paths work
mkdir ignorel
mkdir ignore/figs
mkdir ignore/libraries
mkdir ignore/outputs
mkdir ignore/_targets
mkdir ignore/_targets/controlled
mkdir ignore/_targets/naturalistic
mkdir ignore/stimlists
mkdir ignore/stimuli
mkdir ignore/data

# Pull the matlab folders down from github into ignore/libraries
# If you already have any of these downloaded (say, SPM12)
# You can skip the git clone step and symlink your existing SPM12 folder with ln -s
# so that it shows up at ignore/libraries/spm12 instead
MATLAB_ADDON_PATH=${REPO_PATH}/ignore/libraries
cd $MATLAB_ADDON_PATH
git clone https://github.com/spm/spm12.git
git clone https://github.com/canlab/CanlabCore.git
git clone https://github.com/canlab/Neuroimaging_Pattern_Masks.git
curl -O http://www.animaclock.com/harel/share/gbvs.zip # the source posted by the original author is not on GitHub
unzip gbvs
# extra stuff to compile gbvs
cd gbvs
matlab -nodisplay -r 'gbvs_install; exit'

# recreate conda env
conda env create -f environment.yml
# R package env is handled within R, so not in this setup script

# 2026-01-09 NOT YET AVAILABLE: fMRI DATA NOT ACCESSIBLE YET SORRY!
