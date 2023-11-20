# various helper functions for converting psychopy data to SPM-able events ----

parse_events_naturalistic <- function (file_path) {
  raw <- read_csv(file_path)
  # CHANGE THIS WHEN THE LOG FILES EXPORT TRIGGER TIME PROPERLY -cry-
  t_start <- raw$cross_iti.started[1] - 8
  
  onsets <- raw %>% 
    # Drop the last row of ISI before the end-of-run screen
    filter(!is.na(filename)) %>% 
    select(filename, 
           onset_video = video.started, # Video onset
           offset_video = cross_isi.started, # Changing video offset, backed out from subsequent fixation onset (idiotic but it'll work)
           onset_ratings = mood_pleasantness_slider.started # Rating screens onset
           ) %>% 
    mutate(duration_video = offset_video - onset_video,
           duration_ratings = 12) %>% 
    select(-offset_video) %>% 
    pivot_longer(cols = c(starts_with("onset"), starts_with("duration")),
                 names_to = c(".value", "condition"),
                 names_sep = "_") %>% 
    # Each video gets its own condition! Single trial analysis. Bunkus
    # Ratings all go into the same condition though
    mutate(onset = onset - t_start,
           condition = if_else(condition == "video",
                               str_sub(filename, end = -5L),
                               condition)) %>% 
    select(-filename)
  
  return (onsets)
}

parse_events_nback <- function (file_path) {
  raw <- read_csv(file_path)
  
  # CHANGE THIS WHEN THE LOG FILES EXPORT TRIGGER TIME PROPERLY -cry-
  t_start <- raw$instruct_text.started[2] - 8
  
  onsets <- raw %>% 
    # Drop the wait-for-trigger screen and last row of ISI before the end-of-run screen
    filter(!is.na(filename)) %>% 
    select(run_num = run,
           threatening,
           miniblock_num,
           animal_type,
           hemifield,
           direction, 
           onset_video = video.started, # Video onset
           onset_ratings = mood_pleasantness_slider.started # Rating screens onset
    ) %>% 
    # temporary reconstruction:
    # for odd run nums, location is on the odd miniblocks
    # for even run nums, location is on the even miniblocks
    mutate(attended = if_else(run_num %% 2 == miniblock_num %% 2, "hemifield", "animal")) %>% 
    unite("condition", c(attended, animal_type, hemifield, direction), sep = "_", na.rm = TRUE) %>% 
    mutate(onset = coalesce(onset_video, onset_ratings),
           onset = onset - t_start,
           # if it has ONLY attended, it's a ratings trial
           condition = fct_recode(condition, ratings = "animal", ratings = "hemifield"),
           duration = if_else(condition == "ratings",
                              12,
                              2)) %>% 
    select(condition, onset, duration)
  
  # right now, just for the pilot data, I'm delaying this TODO:
  # get the conditions to go by attended and threatening so that we can collapse across runs later
  
  return (onsets)
  
}

format_events_matlab <- function (onsets, raw_path) {
  # expects long by trial with 3 columns: condition, onset, duration
  onsets %<>%
    chop(cols = c(onset, duration)) %>% 
    mutate(duration = round(map_dbl(duration, unique), digits = 3),
           onset_char = map_chr(onset, \(x) x %>% 
                                  round(digits = 3) %>% 
                                  rvec_to_matlab(row = TRUE) %>% 
                                  str_remove(";"))) %>% 
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
  
  # The easiest/laziest way to retain the naming structure hehe
  out_mat_path <- raw_path %>% 
    # do NOT write this in the raw folder anymore
    str_remove("/raw") %>% 
    # replace the file suffixes
    str_remove("raw.csv") %>% 
    paste0("onsets.mat")
  
  matlab_commands <- c(
    rvec_to_matlabcell(onsets$condition, matname = "names", transpose = TRUE),
    rvec_to_matlabcell(onsets$onset_char, matname = "onsets", transpose = TRUE) %>% 
      # extra formatting idiocy bc these rvec_to_matlab* functions weren't designed to nest
      str_remove_all("'"),
    rvec_to_matlabcell(onsets$duration, matname = "durations", transpose = TRUE) %>% 
      str_remove_all("'"),
    call_function("save",
                  args = c(out_mat_path,
                           "names", "onsets", "durations") %>% 
                    map_chr(wrap_single_quotes))
  )
  
  # we do need out_path to come out for targets tracking
  # even though it's already in the save part of matlab_commands
  return (list(out_path = out_mat_path, matlab_commands = matlab_commands))
}
