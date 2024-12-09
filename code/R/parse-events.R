# various helper functions for converting psychopy data to SPM-able events ----

## actually calculating event timing for onsets ----

parse_events_naturalistic <- function (file_path, conditions_base, tr_duration, n_trs) {
  raw <- read_csv(file_path)
  # in this PsychoPy task, the routines are set up where ITI occurs "before" each video
  # so the first ITI is triggered by the scanner and covers the 8 second discard period
  t_start <- raw$cross_iti.started[2] + (8 %/% tr_duration)*tr_duration
  
  onsets <- raw %>% 
    # As of sub-0001, cross_iti appears "before" its video
    # so video offset must be backed out from the fixation in the next row
    # so for the edge case of the final video, grab cross_iti_final instead
    mutate(offset_video = lead(coalesce(cross_iti.started, cross_iti_final.started))) %>% 
    select(video_id, 
           onset_video = video.started,
           offset_video
           ) %>% 
    mutate(duration_video = offset_video - onset_video,
           duration_ratings = 12) %>% 
    # Now drop the last row of ISI before the end-of-run screen
    filter(!is.na(video_id)) %>% 
    select(-offset_video) %>% 
    pivot_longer(cols = c(starts_with("onset"), starts_with("duration")),
                 names_to = c(".value", "condition"),
                 names_sep = "_") %>% 
    # Each video gets its own condition! Single trial analysis. Bunkus
    # Ratings all go into the same condition though
    mutate(onset = onset - t_start,
           condition = if_else(condition == "video",
                               str_sub(video_id, end = -5L),
                               condition)) %>% 
    # we need to get information from the main stimlist
    # in order to patch it into the condition names
    # to send it into the SPM.mat so that contrasts can be set in matlab
    left_join(conditions_base %>% select(video, animal_type, has_loom), 
              by = c("video_id" = "video")) %>% 
    mutate(condition = glue("{animal_type}_loom{has_loom}_{condition}")) %>% 
    # Just in case there are any onsets that are after the scan actually ended
    # remove them because they will end up creating zero-regressors later which we don't want
    filter(onset + duration < tr_duration * n_trs) %>% 
    select(-video_id)
  
  return (onsets)
}

parse_events_nback <- function (file_path, tr_duration, n_trs) {
  raw <- read_csv(file_path)
  
  # IMPORTANT: At the SPM modeling stage, the first discarded volumes (approx 8 s)
  # are TOTALLY IGNORED FROM THE TIMESERIES
  # so time 0 should ideally be immediately after the PsychoPy 8 s wait
  # In the event that the TR length does not divide evenly into 8,
  # we still need time 0 to line up with the first kept TR,
  # which may be slightly less than 8 s
  # so calculate the number of skipped TRs and get the wait offset from that
  t_start <- raw$wait_screen.started[2] + (8 %/% tr_duration)*tr_duration
  
  onsets <- raw %>% 
    # Drop the wait-for-trigger screen and last row of ISI before the end-of-run screen
    filter(!is.na(miniblock_file)) %>% 
    select(run_num = run,
           threatening,
           attended,
           miniblock_num,
           animal_type,
           hemifield,
           direction, 
           onset_video = video_resp.started, # Video onset. Use the response onset bc video.started doesn't take the 1 s pause into account
           onset_ratings = mood_pleasantness_slider.started # Rating screens onset
    ) %>% 
    group_by(miniblock_num) %>% 
    mutate(attended = attended[1]) %>% 
    ungroup() %>% 
    mutate(attended = case_match(attended, "LOCATION" ~ "hemifield", "ANIMAL" ~ "animal")) %>% 
    unite("condition", c(attended, animal_type, hemifield, direction), sep = "_", na.rm = TRUE) %>% 
    mutate(onset = coalesce(onset_video, onset_ratings),
           onset = onset - t_start,
           # if it has ONLY attended, it's a ratings trial
           condition = fct_recode(condition, ratings = "animal", ratings = "hemifield"),
           duration = if_else(condition == "ratings",
                              12,
                              1.5)) %>% 
    filter(onset + duration < tr_duration * n_trs) %>% 
    select(condition, onset, duration)
  
  # right now, just for the pilot data, I'm delaying this TODO:
  # get the conditions to go by attended and threatening so that we can collapse across runs later
  
  return (onsets)
  
}

realign_onset_end_spike <- function (onsets, add_realign = 0) {
  # only realigns the task stimuli, not the ratings
  # optionally realigns with an additional offset (e.g. to estimate FIR before the offset)
  out <- onsets %>% 
    mutate(offset = onset + duration + add_realign) %>% 
    mutate(onset2 = if_else(condition == "ratings",
                            onset,
                            offset),
           duration2 = if_else(condition == "ratings",
                               duration,
                               0)) %>% 
    select(condition, onset = onset2, duration = duration2)
  
  return (out)
}

# for the FIR analysis, we will just look at looming vs non looming
# so that we can try to power up the basis functions
# so that means that the conditions need to be relabeled at the onsets/level 1 phase
# and not condensed at the contrast phase
relabel_onset_looming_naturalistic <- function (onsets, conditions_base) {
 # only naturalistic needs conditions_base
  # bc nback has this information in the video names
    conditions <- conditions_base %>% 
      select(video, has_loom) %>% 
      mutate(video = str_remove(video, ".mp4"))
    
    out <- onsets %>% 
      left_join(conditions,
                by = c("condition" = "video")) %>% 
      mutate(condition2 = fct_recode(as.character(has_loom),
                                         "looming" = "1",
                                         "receding" = "0"),
             condition2 = fct_relevel(condition2, "looming", "receding")) %>% 
      select(condition = condition2, onset, duration)
    
    return (out)
}

relabel_onset_looming_nback <- function (onsets) {
  onsets %>% 
    mutate(condition = as.character(condition), 
           condition2 = case_when(endsWith(condition, "looming") ~ "looming", 
                                  endsWith(condition, "receding") ~ "receding", 
                                  TRUE ~ "ratings"),
           condition2 = fct_relevel(condition2, "looming", "receding")) %>% 
    select(condition = condition2, onset, duration)
}

format_events_matlab <- function (onsets, raw_path, onset_type) {
  # expects long by trial with 3 columns: condition, onset, duration
  onsets %<>%
    mutate(across(where(is.numeric), \(x) round(x, digits = 3))) %>% 
    chop(cols = onset) %>% 
    # So that across runs, the conditions get output in the same order
    arrange(condition)
  
  # Then it has to be in a Matlab-able format... S A D !!!
  # Per the SPM 12 docs, a single .mat file should contain 3 cell arrays of length == n_conditions
  # Cell array 1: names: each cell contains a string
  # Cell array 2: onsets: each cell contains a vector
  # Cell array 3: durations: each cell contains a vector or a single number if same duration for all trials
  
  # I think it eventually has to be a cell (or struct, that would let you name the conds) array
  # Where each L1 cell corresponds to a condition
  # and then within each cell there's a column for onsets and a column for durations?
  
  out_suffix <- glue("model-{onset_type}_onsets.mat")
  # The easiest/laziest way to retain the naming structure hehe
  out_mat_path <- raw_path %>% 
    # do NOT write this in the raw folder anymore
    str_remove("/raw") %>% 
    # replace the file suffixes
    # use str_sub now to strip off the psychopy date-time suffix
    # which should be formatted with a fixed nchar
    str_sub(end = -28L) %>% 
    str_c(out_suffix)
  
  matlab_commands <- c(
    rvec_to_matlabcell(onsets$condition, sep = ",", matname = "names"),
    rvec_to_matlabcell(onsets$onset, sep = ",", matname = "onsets"),
    rvec_to_matlabcell(onsets$duration, sep = ",", matname = "durations"),
    call_function("save",
                  args = c(out_mat_path,
                           "names", "onsets", "durations") %>% 
                    map_chr(wrap_single_quotes))
  )
  
  # we do need out_path to come out for targets tracking
  # even though it's already in the save part of matlab_commands
  return (list(out_path = out_mat_path, matlab_commands = matlab_commands))
}

## patching confound cols in to get the corresponding beta.nii indices ----

# firm coding the number of confounds here 
# but you can set it by pulling one sample confound file and counting the number of columns
get_beta_indices <- function (events, n_confounds = 24) {
  events %>% 
    mutate(condition = as.character(condition), 
           beta_num = 1:n() + n_confounds*(run-1), 
           beta_name = sprintf("beta_%04d.nii", beta_num)) %>% 
    select(condition, beta_name)
}