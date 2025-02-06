# various helper functions for converting fmriprepped confounds to SPM-able ----

## get path to fmriprep confounds file from subject/run/task ----

get_raw_confounds <- function (subject, task, run) {
  inject(here::here(!!!path_here_derivatives, subject, "func",
                    paste(subject, task, run, "desc-confounds_timeseries.tsv", sep = "_")))
}

## targets-friendly: pull just the realignment motion parameter confounds ----
get_confounds_motion <- function (file_path, tr_duration, disdaq_duration) {
  confounds <- read_confounds(file_path = file_path, 
                              tr_duration = tr_duration, 
                              disdaq_duration = disdaq_duration) %>% 
    # I am pretty sure from inspecting Phil's version of the trimmed confounds
    # that we are including: starts with trans, rot
    # 2024-07-22: Phil said to stop including the spike regressors for motion_outlier
    select(starts_with("trans"), starts_with("rot"))
  # just leave in the outlier cols that are all 0s after filtering off the first volumes
  # so that the columns are standardized
  # SPM will ignore inestimable design matrix columns! nice.
  
  out_path <- file_path %>% 
    # even though matlab load READS comma separated files as arrays just fine
    # when the extension is .csv
    # spm_run_fmri_spec throws a fit if this isn't .txt ending
    str_replace("timeseries.tsv", "trimmed.txt")
  
  # it would seem to be that matlab does not wish these to have headers... sad
  confounds %>% 
    write_csv(file = out_path, col_names = FALSE)
  
  return (out_path)
}

## targets-friendly: estimate dvars outlierness for a single run ----
get_noise_measures <- function (file_path, tr_duration, disdaq_duration) {
  confounds <- read_confounds(file_path = file_path, 
                              tr_duration = tr_duration, 
                              disdaq_duration = disdaq_duration) %>% 
    select(contains("global_signal"), contains("dvars"), "framewise_displacement")
  
  return (confounds)
}

## helper: read in fmriprep confounds file as tibble, discard the first volumes ----
read_confounds <- function (file_path, tr_duration, disdaq_duration) {
  confounds <- read_tsv(file_path, na = "n/a")
  n_trs <- nrow(confounds)
  # nb: int division makes it work even if the TR does not divide evenly into 8
  # int div is floor div so this may start at a TR slightly before 8 s has elapsed
  start_tr <- (disdaq_duration %/% tr_duration) + 1
  
  out <- confounds %>% 
    # get the TR programmatically pls!!!
    # remove the first 8 seconds of data
    slice(start_tr:n())
  
  return (out)
}
