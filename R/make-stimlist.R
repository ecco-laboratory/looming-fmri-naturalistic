# various helper functions for making stimlist-related targets ----

nback_block_threatenings <- c("dog", "frog", "spider", "above", "below")

make_stims_nback <- function (base_stim_block,
                              threatening = NULL,
                              predictable = NULL,
                              n_stimuli_per_condition = 6,
                              p_threatening = 2/3) {

  # operates on one block's worth of stimuli
  
  if (!is.null(threatening) & !is.null(predictable)) {
    out <- base_stim_block %>% 
      # this is a temp helper column for generating the stim-specific loom/recedes
      mutate(threat_status = case_when(
        (predictable == "animal_type" & animal_type == threatening) | (predictable == "hemifield" & hemifield == threatening) ~ "threatening",
        (predictable == "animal_type" & animal_type != threatening) | (predictable == "hemifield" & hemifield != threatening) ~ "safe",
        TRUE ~ NA_character_),
        direction = map(threat_status, \(x) {
          if (x == "threatening") {
            directions <- rep(c("looming", "receding"), times = c(n_stimuli_per_condition*p_threatening,
                                                                  n_stimuli_per_condition*(1-p_threatening)))
          } else {
            directions <- rep(c("looming", "receding"), times = c(n_stimuli_per_condition*(1-p_threatening),
                                                                  n_stimuli_per_condition*p_threatening))
          }
          return (directions)
        })
        ) %>%
      select(-threat_status)
    
  } else if (is.null(threatening) & is.null(predictable)) {
    out <- base_stim_block %>% 
      mutate(
        # doesn't actually depend on the column content so just map anything
        direction = map(animal_type, \(x) {
          rep(c("looming", "receding"), each = n_stimuli_per_condition*0.5)
        })
      )
      
  } else {
    stop("Threatening and predictable need to BOTH be set or BOTH unset!")
  }
  
  out %<>%
    unchop(direction) %>% 
    # trial ordering will be computed by the canlabtools matlab script
    # this select really just reorders the cols for convenience
    select(attended, animal_type, hemifield, direction)
  
  return (out)
}

make_stimlist_nback <- function (out_path_glue,
                                 trial.order_path,
                                 metadata,
                                 threatening = NULL,
                                 predictable = NULL
                                 ) {
  
  onsets <- read_csv(trial.order_path, col_names = "condition_num") %>% 
    mutate(condition_num = as.integer(condition_num),
           tr_num = 1:n())
  
  # filename is going to be reserved in psychopy for writing out the data
  # so don't overload it, just to be safe
  metadata %<>%
    rename(video_id = filename)
  
  if (!is.null(threatening) & !is.null(predictable)) {
    # NOT putting "attended" in here
    # it will be shown by PsychoPy at the beginning of each block
    # but doesn't otherwise change the shown stimuli
    # the threatening category varies between scan run
    # append it into the files here to make analysis easier later
    onsets %<>%
      mutate(predictable = predictable,
             threatening = threatening)
  } else if (is.null(threatening) & is.null(predictable)) {
  } else {
    stop("Threatening and predictable need to BOTH be set or BOTH unset!")
  }
  
  tr_max <- nrow(onsets)
  
  onsets %<>%
    mutate(miniblock_num = (tr_num-1) %/% (tr_max/nback_blocks_per_run) + 1) %>% 
    # do this to label the rating rest TRs with another condition_num
    # so that they don't get counted by the ITI duration calculator below
    group_by(miniblock_num) %>% 
    mutate(miniblock_tr_num_rev = n():1,
           condition_num = if_else(miniblock_tr_num_rev <= 12,
                                   13L,
                                   condition_num)) %>% 
    select(-miniblock_tr_num_rev) %>% 
    ungroup() %>% 
    # get it to one row per fixation period
    filter(condition_num + lag(condition_num) > 0) %>% 
    # bind on stimulus filenames for each condition
    nest(trials = -condition_num) %>% 
    mutate(trials = map2(trials, condition_num, \(x, y) {
      # nothing 2 bind if it's the ITI condition
      if (y %in% c(0, 13)) return (x)
      metadata_to_bind <- metadata %>% 
        filter(condition_num == y) %>% 
        slice_sample(n = nrow(x), replace = TRUE) %>% 
        # must drop condition_num before binding
        # or else the column names will conflict
        # can't join it because there are no row-unique cols
        select(-condition_num)
      
      return (bind_cols(x, metadata_to_bind))
    })) %>% 
    unnest(trials) %>% 
    # nesting by trials changes the row order, get it back in chronological
    arrange(tr_num) %>% 
    # pull up the subsequent ITI duration
    # into a column next to its preceding trial
    mutate(iti = if_else(condition_num > 0 & lead(condition_num == 0),
                         (lead(tr_num, 2) - lead(tr_num)) * 2,
                         NA_integer_)) %>% 
    # drop the extraneous ITI and rating rest rows now
    filter(condition_num > 0, condition_num < 13) %>%
    # patch in 0 otherwise (will work for consecutive stims and final stims)
    mutate(iti = coalesce(iti, 0),
           trial_num = 1:n())
  # TODO: Deal with leading ITI. Like wuzzthedeal
  # Get 1 csv per miniblock. Seems necessary :(
  
  onsets %>%
    nest(trials = -miniblock_num) %>% 
    mutate(out_file = glue::glue(out_path_glue)) %$% 
    # to write the files out
    # walk2 should invisibly return out_file for targets
    walk2(out_file, trials,
          \(x, y) {write_csv(y,
                             file = x)})
}
