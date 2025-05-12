# various helper functions for converting psychopy data to SPM-able events ----

## get path to events file from subject/run/task ----

# run should be in run-%02d bids format
get_raw_events <- function (subject, task, run) {
  file <- list.files(here::here("ignore", "data", "beh", subject, "raw"),
                      pattern = paste(subject, task, run, sep = "_"),
                      full.names = TRUE)
  # need it to explicitly error if the raw events haven't been uploaded yet
  # or if multiple matching files are found (e.g. from re-starting a psychopy task)
  stopifnot(length(file) == 1)
  file
}

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
  
  if (onset_type == "endspike") onsets %<>% realign_onset_end_spike()
  
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
  
  # go ahead and run the matlab code right from here
  # assume matlab_path is a global variable that will be instantiated in the targets script
  out <- run_matlab_target(matlab_commands, out_mat_path, matlab_path)
  
  return (out)
}

## constructing encoding model activation timecourses from onset orders ----

# this one just transforms an onsets-formatted dataframe into one with 1 row per TR
# and a column saying what stimulus (or fixation) is on screen
make_condition_timecourse <- function (onsets,
                                       out_path,
                                       tr_duration, 
                                       n_trs) {
  stims <- tibble()
  
  for (i in 0:(n_trs-1)) {
    this_time <- i * tr_duration
    # the most recent stimulus to have appeared on screen (might have already gone away)
    # bind 0 onto the indices to handle the case of TRs before the first stim has appeared
    idx_last_stim <- max(c(0L, which(this_time - onsets$onset >= 0)))
    # this labels the entire TR as whatever's on screen at the beginning of the TR
    if (idx_last_stim > 0) {
      stim_still_onscreen <- this_time - onsets$onset[idx_last_stim] <= onsets$duration[idx_last_stim]
    } else {
      stim_still_onscreen <- FALSE
    }
    
    if (stim_still_onscreen) {
      this_stim <- onsets %>% 
        slice(idx_last_stim) %>% 
        select(condition) %>% 
        mutate(tr_num = i+1)
    } else {
      # if not, it's fixation
      this_stim <- tibble(tr_num = i+1,
                          condition = "fixation")
    }
      
    stims <- stims %>% 
      bind_rows(this_stim)
  }
  
  write_csv(stims,
            file = out_path)
  
  return (out_path)
}

# takes single-subject timecourse (result of make_condition_timecourse but row-concatenated across runs)
# and relabels the video_id column not to correspond to the actual presentation,
# but an endspike window anchored around the END of each video.
# using R convention, 1 corresponds to the first TR of each video and -1 corresponds to the final TR,
# so you can count either forwards or backwards
# window length is inclusive of start and end
# returns with an appropriate lag-windowed version of the timecourse
# where TRs not belonging to a stimulus are marked as NA
relabel_timecourse_endspike <- function (timecourse,
                                         window_start = -2L,
                                         window_length = 10L) {
  
  timecourse$video_id_lagged <- NA
  
  for (i in 1:(nrow(timecourse)-1)) {
    # if window_start is negative, count from the END of the video BACKWARDS
    if (sign(window_start) == -1) {
    # if a final TR is identified, label a 10-TR (4.92 second) window from the last 2 TRs within the video to the 8 TRs after
    # to get end-of-video related BOLD variation??? peut-etre?
    if (timecourse$video_id[i] != "fixation" & timecourse$video_id[i+1] == "fixation") {
      this_window_start <- i+1+window_start
      this_window_end <- this_window_start+window_length-1
      timecourse$video_id_lagged[this_window_start:this_window_end] <- timecourse$video_id[i]
      }
    } else {
      # otherwise if positive, count from the BEGINNING of the video FORWARDS
      # identify initial TRs only
      if (timecourse$video_id[i] != "fixation" & timecourse$video_id[i-1] == "fixation") {
        this_window_start <- i-1+window_start
        this_window_end <- this_window_start+window_length-1
        timecourse$video_id_lagged[this_window_start:this_window_end] <- timecourse$video_id[i]
      }
    }
  }
  
  out <- timecourse %>% 
    select(-video_id) %>% 
    rename(video_id = video_id_lagged)
  
  return (out)
}

# similar to above, but instead of returning fixed-length windows anchored at different points in the stimulus video
# returns the full original boxcar but with an added tail of "stimulus"-labeled TRs at the end
# to account for hemodynamic lag relative to stimulus presentation
# tail length of 0 is equivalent to the original timecourse
# it doesn't handle negative tails right now so don't try it
relabel_timecourse_boxcartail <- function (timecourse,
                                           tail_length = 10L) {
  timecourse$video_id_lagged <- NA
  
  for (i in 1:(nrow(timecourse)-1)) {
    # if it's any video TR, label it
    if (timecourse$video_id[i] != "fixation") {
      timecourse$video_id_lagged[i] <- timecourse$video_id[i]
      # if it's the final TR, calculate and label the tail as well
      if(timecourse$video_id[i+1] == "fixation") {
        this_window_end <- i+tail_length
        timecourse$video_id_lagged[i:this_window_end] <- timecourse$video_id[i]
      }
    }
  }
  
  out <- timecourse %>% 
    select(-video_id) %>% 
    rename(video_id = video_id_lagged)
  
  return (out)
}

# this one assembles activation timecourses and then outputs to header-less tabular data for matlab
make_encoding_timecourse_matlab <- function (onsets, 
                                             path_stim_activations,
                                             out_path,
                                             tr_duration,
                                             run_duration,
                                             fixation_activation = NULL, 
                                             stim_fps = 10L) {
  
  # these come in from python so the unit numbers in the col names start indexing at 0
  stim_activations <- read_csv(path_stim_activations, name_repair = "unique") %>% 
    select(-frame) %>% 
    nest(activations = -video)
  
  onsets %<>%
    mutate(video = paste0(str_sub(condition, start = -7L), ".mp4"))
  
  if (!is.null(fixation_activation)) {
    # use the 'veridical' activations from the fixation video through FlyNet or wherever
    # 2024-12-16: NOT RECOMMENDED BY MONICA bc the baseline activation for zero-motion stimuli is kind of high
    # and sure, maybe somewhere like SC is doing defaulty stuff where task activation actually deactivates it
    # but maybe it isn't??
  } else {
    # and then set the activations during fixation to the mean of each unit's activation to all the stims in this run
    fixation_activation <- stim_activations %>% 
      semi_join(onsets, by = "video") %>% 
      unnest(activations) %>% 
      summarize(across(-video, mean)) %>% 
      mutate(video = "fixation")
  }
  
  # fencepost first condition! 
  # there may be some tiny amount of included fixation before the very first onset
  # if it's more than 1 frame's worth of time we should include it
  if (floor(onsets$onset[1] * stim_fps) > 0) {
    timecourse <- rep(list(fixation_activation), times = floor(onsets$onset[1] * stim_fps)) %>% 
      bind_rows() 
  } else {
    timecourse <- tibble()
  }

  for (i in 1:nrow(onsets)) {
    these_activations <- stim_activations %>% 
      filter(video == onsets$video[i]) %>% 
      # unnest to keep the video column. in case it's useful later!
      unnest(activations)
    
    # if the current stim isn't in the dictionary of all stim activations, something is wrong
    stopifnot(nrow(these_activations) > 0)
    
    timecourse <- timecourse %>% 
      bind_rows(these_activations)

    # calculate the time from the current end of the timecourse to the next onset
    # and generate/append the requisite amount of fixation
    # if it's the very last video, calculate the time to the end of run
    if (i < nrow(onsets)) next_onset <- onsets$onset[i+1] else next_onset <- run_duration
    
    
    this_iti <- next_onset - nrow(timecourse)/stim_fps
    this_fixation_activation <- rep(list(fixation_activation), 
                                    times = round(this_iti * stim_fps)) %>% 
      bind_rows()
    
    timecourse <- timecourse %>% 
      bind_rows(this_fixation_activation)
  }
  
  # linear interpolate every column into TR frequency
  out <- timecourse %>% 
    select(-video) %>% 
    as.list() %>% 
    # approx only works colwise so cheese the map-ing
    map(\(x) approx(1:length(x)/stim_fps, 
                    x, 
                    xout = seq(from = tr_duration, 
                               to = run_duration, 
                               by = tr_duration)) %>% 
          pluck("y")) %>% 
    as_tibble()
  
  write_csv(out,
            file = out_path,
            # because they're going into evil MATLAB
            col_names = FALSE)
  
  return (out_path)
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