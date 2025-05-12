# designed to run on a _single_ parcel-beta
make_rdms_from_beta <- function (path_parcel_beta) {
  
  betas <- read_csv(path_parcel_beta)
  # in case the ROI has 0 voxels? which happens apparently?
  if (nrow(betas) == 0) return (NULL)
  
  roi <- str_sub(basename(path_parcel_beta), start = 7L, end = -5L)
    
  betas %<>%
    # VarX and Row are the names that come in from matlab writetable()
    pivot_longer(cols = starts_with("Var")) %>% 
    pivot_wider(names_from = Row) %>% 
    select(-name) %>%
    cor() %>% 
    as_tibble(rownames = "condition_row") %>% 
    pivot_longer(cols = -condition_row,
                 names_to = "condition_col",
                 values_to = "correlation") %>% 
    # don't forget to put the ROI name into the dataframe for row-binding later!
    mutate(roi = roi)
  
  return (betas)
}

halve_tidy_rdm <- function (rdm_full, row_col, col_col) {
  row_col <- enquo(row_col)
  col_col <- enquo(col_col)
  
  rdm_halved <- rdm_full %>% 
    mutate(sorter = map2_chr({{row_col}}, {{col_col}}, \(x, y) paste(sort(c(x, y)), collapse = " "))) %>% 
    distinct(sorter, .keep_all = TRUE) %>% 
    select(-sorter)
  
  return (rdm_halved)
}
