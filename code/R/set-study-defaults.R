# GET AND SET STUDY DEFAULTS ----
# THIS SCRIPT GETS AND SETS STUDY-WIDE CONSTANTS AND WRITES THEM TO A JSON
# FOR USE BY OTHER SCRIPTS, including in other languages
# this script is not tar_source()-d but regular source()-d bc it's used for target construction
# but not within any targets
# hence we do need to load libraries here
require(RNifti)
require(tidyverse)
require(rlang)

## Set paths to fMRI data ----
path_here_fmri <- c("ignore", "data", "fmri")
path_here_derivatives <- c(path_here_fmri, "derivatives", "fmriprep-23.1.4")

## Pull scan info stored in nifti headers ----
# Assume that sub-0001 has the standard study scans!
print("loading controlled header...")
hdr_controlled <- inject(here::here(!!!path_here_fmri, "sub-0001", "func", "sub-0001_task-controlled_run-01_bold.nii.gz")) %>% 
  readNifti() %>% 
  niftiHeader()

n_trs_controlled <- hdr_controlled$dim[5]
  
tr_duration <- round(hdr_controlled$pixdim[5], digits = 3) # in seconds

print("loading naturalistic header...")
hdr_naturalistic <- inject(here::here(!!!path_here_fmri, "sub-0001", "func", "sub-0001_task-naturalistic_run-01_bold.nii.gz")) %>% 
  readNifti() %>% 
  niftiHeader()

n_trs_naturalistic <- hdr_naturalistic$dim[5]

# both tasks should have the same TR
stopifnot(round(hdr_naturalistic$pixdim[5], digits = 3) == tr_duration)

## pull info from scans.tsv ----

n_task_runs <- inject(here::here(!!!path_here_fmri, "sub-0001", "sub-0001_scans.tsv")) %>% 
  read_tsv() %>% 
  separate_wider_delim(cols = filename, delim = "/", names = c("folder", "filename")) %>% 
  filter(folder == "func") %>% 
  separate_wider_delim(cols = filename, delim = "_", names = c(NA, "task", "run", NA)) %>% 
  count(task) %>% 
  rename(n_runs = n)

## calculate info related to the disdaq wait period at the beginning of each run ----
disdaq_duration <- 8 # also in seconds! CONSTANT PER PHIL! Set by experimenters (us)
# use floor division to calculate the number of TRs discarded 
# to get as close to the disdaq duration as possible without going over
# the scanner tasks wait for the full disdaq duration 
# so the first kept TR should start slightly before the wait period ends
n_trs_kept_controlled <- n_trs_controlled - (disdaq_duration %/% tr_duration)
# for pilots sub-9902-9904, this was 1407
# but before sub-0001, removed in-scanner ratings so the run is now shorter
n_trs_kept_naturalistic <- n_trs_naturalistic - (disdaq_duration %/% tr_duration)

study_defaults <- list(path_fmri = inject(here::here(!!!path_here_fmri)),
     path_fmri_derivatives = inject(here::here(!!!path_here_derivatives)))

task_defaults <- n_task_runs %>% 
  mutate(n_trs = if_else(endsWith(task, "controlled"), n_trs_controlled, n_trs_naturalistic),
         n_trs_kept = if_else(endsWith(task, "controlled"), n_trs_kept_controlled, n_trs_kept_naturalistic),
         tr_duration = tr_duration,
         disdaq_duration = disdaq_duration)

study_defaults %>% 
  jsonlite::write_json(path = here::here("study_defaults.json"))

task_defaults %>% 
  jsonlite::write_json(path = here::here("task_defaults.json"))
