# Human superior colliculus pathways represent the form and motion of looming objects

The superior colliculus is well known for its roles in visual orienting, oculomotor control, attention, and defensive behavior across species. Recently, [we predicted and found](https://www.cell.com/iscience/fulltext/S2589-0042(24)01108-8) that representations from a shallow convolutional neural network could predict defensive blinking to looming objects in infants and superior colliculus responses to optical expansion in adults. These findings suggest that the superior colliculus may coordinate defensive responses to looming in humans. 

In this project, we tested whether the human superior colliculus functions in isolation during looming threat perception, or if it covaries with cortical networks involved in visual salience and object recognition. We used computational models of [looming detection](https://elifesciences.org/articles/72067), [visual saliency](https://ieeexplore.ieee.org/document/730558/), and [object recognition](https://papers.nips.cc/paper_files/paper/2012/hash/c399862d3b9d6b76c8436e924a68c45b-Abstract.html) to predict patterns of superior colliculus BOLD response acquired as participants viewed a series of naturalistic videos or performed a working memory task, and to examine their functional connectivity with the rest of the brain. A detailed description of this work is forthcoming. This repository provides the source code used to conduct these analyses.

<p align="center">
<img src="https://github.com/ecco-laboratory/looming-fmri-naturalistic/blob/main/GraphicalAbstract.png" width="400">
</p>

Please post any questions about the analysis repository as repo issues and we'll get to you as soon as we can. Enjoy!

## SHORT VERSION

**Follow these steps to install dependencies and run code.**

1. Clone the repo
1. Install the following:
    1. MATLAB, ideally no older than R2024a
    1. datalad
    1. conda
    1. R 4.4.1 and RStudio
1. Edit the following:
    1. environment.yml: Set the env name and prefix to match the local path to your clone folder
    1. R/_targets_retinotopy.R: edit the `matlab_path` variable to point to your own Matlab install path
1. Run the following:
    1. In a terminal: `code/shell/init_project.sh`
    1. In RStudio (opened using the .Rproj file): `renv::restore()`, follow all prompts as necessary
    1. In the same RStudio, `source('run.R')` to **run the whole analysis pipeline.**

## LONG VERSION

### Pull code

First, clone this repository to your favorite location on your computer. :3 The setup scripts and analysis pipelines rely almost exclusively on setting things up in subfolders of the cloned project repo.

### Install/edit some dependencies manually

Next, there are a couple things you'll have to set up manually before you can start running my helper scripts.

#### Some other programs & command line tools

The pipeline requires the following programs and command line tools that I can't really auto-download for you. Please download the following if you don't have them already:

- **MATLAB:** The fMRI analyses use SPM and SPM-dependent tools. The analyses were originally written on MATALB R2024a.
- [**conda:**](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) used for Python analyses (the core neural network model is implemented through PyTorch) and maintaining the Python package environment
- **[R](https://cran.r-project.org) (ideally 4.4.1) and [RStudio:](https://posit.co/download/rstudio-desktop/)** used for R analyses (all tabular data analysis, major statistics, and graphs across 3 studies)

#### Manually edit some file paths

I have done everything possible to have the setup scripts and analysis pipeline set up everything to relative path to specific subfolders of the project repo, so that you need to edit as few file paths as possible. However, there were still a couple paths I couldn't totally automate that you will need to manually change before you continue:

1. In `environment.yml`, set the env name and prefix to match the local path to your clone folder, or wherever you want to create the conda environment.
1. In `.Rprofile`, where global variables are set for the R targets pipeline:
    1. Edit the `conda_path` variable to point _directly to the binary_ associated with the environment you cloned from our `environment.yml`. (It should point to the `/bin` inside the conda env folder.)
    1. Edit the `matlab_path` variable to point to your own Matlab install path.

### Run first package setup script

I have included a `code/shell/init_project.sh` script that _should_ download the following pipeline dependencies for you. If you want to know what this setup script is doing, read on!

#### Download additional Matlab libraries

- [SPM12](https://github.com/spm/spm12)
- [canlab/CanlabCore](https://github.com/canlab/CanlabCore)
- [canlab/Neuroimaging_Pattern_Masks](https://github.com/canlab/Neuroimaging_Pattern_Masks)
- [GBVS toolbox](http://www.animaclock.com/harel/share/gbvs.php)

These are all distributed as folders of functions that have no installer beyond cloning the repo to your local and then adding the folder of functions to your Matlab path. The `shell/init_project.sh` setup script should `git clone` all of the libraries into the `ignore/libraries` subfolder for you. All of the Matlab scripts point to these local install paths, so they should then find the folders correctly without additional modification!

If you have one or more of these libraries already (most likely SPM12), you can change the setup script to use `ln -s` (on Mac or Linux) to symlink your existing SPM12 folder to appear under `ignore/libraries/spm12/` as well, instead of git cloning a duplicate copy of SPM12.

#### Set up Python packages (with conda)

We have included a conda `environment.yml` file that should automatically install Python 3.9.13 (the version we used) and the Python package dependencies necessary to reproduce our analyses. `code/shell/init_project.sh` calls 

```python
conda env create -f environment.yml
```

for you.

#### Download fMRI data with datalad

**NOT YET AVAILABLE!** As soon as it is, we will update `init_project.sh` to download the fMRI dataset from OpenNeuro for you. Please mind that it will download into a subfolder of this project folder! If you want to download the large dataset to another folder (perchance on another drive more suited to large file storage), you can similarly change the setup script to use `ln -s` symlink the dataset to appear under `ignore/data/fmri`.

### Set up R packages (with renv)

We have included an R project file and `renv` package environment lockfiles that should automatically install the R package dependencies necessary to reproduce our analyses.

We used R 4.4.1--renv will not manage the associated R version for you, only the packages, so we recommend you make sure you have this version installed before recreating our repository. renv does not _require_ you to have the specific R version encoded in the lockfile in order to restore a package environment, but I _strongly recommend_ you use the specific R version, which should allow renv to pull pre-compiled binaries for all package dependencies. You can make it work with a slightly newer R version, but renv may attempt to install certain packages from C++ source, which can cause you headaches.

If you have a different (probably newer) version of R installed, follow [these instructions](https://jacobrprice.github.io/2019/09/19/Installing-multiple-parallel-R-versions.html) to set up multiple side-by-side R versions on your Mac.

Follow [these instructions](https://support.posit.co/hc/en-us/articles/200486138-Changing-R-versions-for-the-RStudio-Desktop-IDE) to choose the active R version. For Mac users, I recommend installing and using the RSwitch menu bar utility.

Once you have R 4.4.1 (or whatever newer R version if you want to gamble with fate):

1. Open the .Rproj file in this repo. Once RStudio is open to the project folder, the renv autoloader should run on startup to install the renv package if you don't have it, and then the renv environment will activate.
1. Call `renv::restore()` to recreate the package environment.
1. Follow prompts as necessary to install all packages in the environment.

### Run pipeline through R

I use the [targets](https://books.ropensci.org/targets/) package to manage dependency tracking for data inputs, processing scripts, and outputs. (For example, the pipeline uses the `osfr` R package to automatically download task stimuli etc. off of OSF for you and save it into a standard subfolder of the project folder.)

I have declared a separate targets pipeline for each study (Study 1: naturalistic, Study 2: controlled) in two respective `code/R/_targets_[THIS_STUDY].R` scripts.

You can run the `run.R` script as follows:

```bash
cd path/to/looming-fmri
Rscript -e 'source("run.R")'
```

`run.R` attempts to call `targets::tar_make()` to execute all of the analysis code for Study 2 and then 1. (I have them written in that order because there are some Study 1 analyses that take Study 2 results as dependencies.) Hopefully it'll just... work!

If you have slurm set up for your computing environment, I have set up the `_targets_*.R` scripts to run via slurm by default, to attempt to parallel the analysis operations. Please read the extra notes I've given in `run.R` if you want to attempt slurm, because I fear some of it is environment-specific.  I've tried to document everything I can but I make no promises about the slurm controller working for you out of the box!

