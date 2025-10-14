# TODO: Make sure your qualtrics API token is reproducible yet safe from prying eyes
# TODO also: set up another target so that targets will think the survey is updated when it is
# This helper function parses this specific questionnaire, so the survey ID is hard-coded within, not an arg
get_splat_stimulus_norms_qualtrics <- function () {
  qualtrics_data <- fetch_survey("SV_6nZOkqUU0MEAnMa",
                                 include_display_order = FALSE)
  
  video_qualtrics_ids <- qualtrics_data %>% 
    select(contains("caption")) %>% 
    sjlabelled::get_label() 
  
  df_video_qualtrics_ids <- tibble(video_num = names(video_qualtrics_ids),
                                   qualtrics_id = video_qualtrics_ids %>% 
                                     str_split_i(pattern = " ", i = 1)) %>% 
    separate_wider_regex(cols = video_num,
                         patterns = c(video_num = "[[:digit:]]+",
                                      "_[[:alpha:]]+",
                                      block_num = "[[:digit:]]+")) %>% 
    mutate(across(ends_with("num"), as.integer))
  
  out <- qualtrics_data %>% 
    filter(Status != "Survey Preview", participantId_check != "test") %>% 
    select(-contains("quit"), -starts_with("Recipient"), -DistributionChannel, -UserLanguage) %>% 
    pivot_longer(cols = c(contains("caption"), 
                          contains("pleasantness"), 
                          contains("arousal"), 
                          contains("fear")), 
                 names_to = c("video_num", ".value", "block_num"), 
                 names_pattern = "([[:digit:]]+)_([[:alpha:]]+)([[:digit:]]+)",
                 names_transform = list(video_num = as.integer,
                                        block_num = as.integer)) %>% 
    left_join(df_video_qualtrics_ids, by = c("video_num", "block_num"))
  
  return (out)
}

join_raw_norms_to_stim_labels <- function (norms_raw, annotations, loom_col, animal_col = "animal_type") {
  out <- inner_join(norms_raw, 
                    annotations,
                    by = "qualtrics_id") %>% 
    select(participantId, 
           video_id, 
           !!loom_col, 
           !!animal_col, 
           caption, 
           rating_pleasantness = pleasantness, 
           rating_arousal = arousal, 
           rating_fear = fear) %>% 
    filter(!is.na(rating_pleasantness), !is.na(rating_arousal), !is.na(rating_fear))
  
  return (out)
}