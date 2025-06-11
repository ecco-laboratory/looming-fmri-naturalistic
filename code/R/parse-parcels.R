get_parcel_conjunctions <- function (path_parcels_1, 
                                     path_parcels_2, 
                                     threshold_p = .05,
                                     p_adjust_method = "BH") {
  tvals1 <- path_parcels_1 %>% 
    get_parcel_tvals_long() %>% 
    label_parcel_pvals_long()
  
  tvals2 <- path_parcels_2 %>% 
    get_parcel_tvals_long() %>% 
    label_parcel_pvals_long()
  
  out <- full_join(tvals1, tvals2, by = "label", suffix = c("_1", "_2")) %>% 
    mutate(tval_conj = pmin(tval_1, tval_2),
           pval_conj = map_dbl(tval_conj, \(x) map_pval_from_tval(x)),
           pval_conj = p.adjust(pval_conj, method = p_adjust_method, n = n()))
  
  return (out)
}

get_parcel_diff_tvals <- function (path_parcels_1, 
                                   path_parcels_2) {
  conn1 <- get_parcel_values_long(path_parcels_1)
  
  conn2 <- get_parcel_values_long(path_parcels_2)
  
  out <- full_join(conn1, conn2, by = c("fold_num", "label"), suffix = c("_1", "_2")) %>% 
    mutate(diff_connectivity = value_1 - value_2) %>% 
    group_by(label) %>% 
    summarize(tval = mean(diff_connectivity)/(sd(diff_connectivity)/sqrt(n())))
  
  return (out)
}

label_parcel_pvals_long <- function (parcels_long, 
                                     threshold_p = .05, 
                                     p_adjust_method = "BH") {
  out <- parcels_long %>% 
    mutate(pval = map_dbl(tval, \(x) map_pval_from_tval(x)),
           pval = p.adjust(pval, method = p_adjust_method, n = n()))
  
  return (out)
}

get_parcel_tvals_long <- function (path_parcels) {
  out <- get_parcel_values_long(path_parcels) %>% 
    group_by(label) %>% 
    summarize(tval = mean(value)/(sd(value)/sqrt(n())))
  return (out)
}

get_parcel_values_long <- function (path_parcels) {
  out <- path_parcels %>% 
    read_csv() %>% 
    # currently subj nums are not saved out from matlab_parcellate_avg 
    # cause we don't SEEM to need them but that could change
    mutate(fold_num = 1:n()) %>% 
    pivot_longer(cols = -fold_num, names_to = "label", values_to = "value")
  
  return (out)
}

relabel_glasser_clt2ggseg <- function (parcels_long) {
  cortex_parcels <- parcels_long %>% 
    filter(startsWith(label, "Ctx")) %>% 
    mutate(hemi2 = str_sub(label, start = -1L),
           region = str_sub(label, start = 5L, end = -3L),
           region = str_replace(region, "_", "-"),
           hemi1 = if_else(hemi2 == "L", "lh", "rh")) %>% 
    # remove = FALSE to keep region (it's useful later!) and then select to drop hemi1, which is lh/rh
    unite(col = "label", hemi1, hemi2, region, sep = "_", remove = FALSE) %>% 
    select(-hemi1) %>% 
    rename(hemi = hemi2)
    
    
    subcortex_parcels <- parcels_long %>% 
      filter(!startsWith(label, "Ctx")) %>% 
      # because Amygdala doesn't have a prefix but Bstem_SC and Thal_* do
      separate_wider_delim(cols = label, delim = "_", names = c(NA, "region"), too_few = "align_end") %>% 
      mutate(label = paste0("bl_B_", region))
    
    return (bind_rows(cortex_parcels, subcortex_parcels))
}
