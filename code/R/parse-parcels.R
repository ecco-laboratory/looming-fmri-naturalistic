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

label_parcel_pvals_long <- function (parcels_long, 
                                     threshold_p = .05, 
                                     p_adjust_method = "BH") {
  out <- parcels_long %>% 
    mutate(pval = map_dbl(tval, \(x) map_pval_from_tval(x)),
           pval = p.adjust(pval, method = p_adjust_method, n = n()))
  
  return (out)
}

get_parcel_tvals_long <- function (path_parcels) {
  out <- path_parcels %>% 
    read_csv() %>% 
    summarize(across(everything(), \(x) mean(x)/(sd(x)/sqrt(length(x))))) %>% 
    pivot_longer(cols = everything(), names_to = "label", values_to = "tval")
  
  return (out)
}

relabel_glasser_clt2ggseg <- function (parcels_long) {
    parcels_long %>% 
        filter(startsWith(label, "Ctx")) %>% 
        mutate(hemi2 = str_sub(label, start = -1L),
        label = str_sub(label, start = 5L, end = -3L),
        label = str_replace(label, "_", "-"),
        hemi1 = if_else(hemi2 == "L", "lh", "rh")) %>% 
        unite(col = "label", hemi1, hemi2, label, sep = "_")
}
