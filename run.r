#!/usr/bin/env Rscript

# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

# Leaving USE_SLURM as FALSE will build each pipeline _in series,_
# one intermediate object at a time. Nothing wrong with this,
# But it will take a while!

# If you can set up this project on a cluster with slurm, 
# you can reproduce the study analyses in parallel!
# IN ADDITION TO CHANGING USE_SLURM to TRUE, YOU WILL NEED TO:
# Go to the respective _targets_[STUDY].R files
# and edit the slurm parameters being passed to the targets/crew job controller
# in the `tar_options_set()` call near the tops of the targets scripts.
# Refer to the crew.cluster docs to see all the args, but you should recognize the slurm args.
# Set the --account and partition flags to match the names on your cluster.
# Special arg: The workers argument of crew_controller_slurm() controls how many jobs targets
# will attempt to launch in parallel at once.

# If you are experienced with targets, you can even change the controller from crew_controller_slurm()
# to a different crew controller to use your desired engine (including local forking R subprocesses!)
# to launch jobs in parallel.

USE_SLURM=FALSE

Sys.setenv(TAR_PROJECT='controlled') # Study 2
# Run everything referenced in Study 1 analyses or directly in the main text
targets::tar_make(
    c(
        metadata_videos_nback,
        plot_norms,
        beh,
        activations.flynet_raw,
        activations.alexnet_raw,
        starts_with('con')&contains('boxcar'),
        summary_bidsmreye.decoding_fixation_hemifield,
        summary_bidsmreye.decoding_fixation_direction,
        summary_bidsmreye.decoding_fixation_obj
    ),
    use_crew=!USE_SLURM
)
Sys.setenv(TAR_PROJECT='naturalistic') # Study 1
# Run everything referenced in the main text or supplement
targets::tar_make(c(rmd_ms_stats, rmd_ms_supp), use_crew=!USE_SLURM)
