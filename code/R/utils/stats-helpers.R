cohens.d <- function (vec1, vec0, na.rm = FALSE) {
  (mean(vec1, na.rm = na.rm) - mean(vec0, na.rm = na.rm)) / sqrt((var(vec1, na.rm = na.rm)+var(vec0, na.rm = na.rm))/2)
}

cohens.d.2 <- function (mean1, sd1, mean0, sd0) {
  (mean1 - mean0) / sqrt((sd1^2+sd0^2)/2)
}

softmax <- function (vec) {
  denom <- sum(exp(vec))
  return (exp(vec) / denom)
}

summary_stats_default <- list(mean = mean,
                              ci95.lower = \(x) quantile(x, .025),
                              ci95.upper = \(x) quantile(x, .975))

map_pval_from_tval <- function (x, two.tailed = TRUE) {
  if (is.nan(x)) {
    return (NaN)
  } else {
    if (two.tailed) {
        return (pnorm(x, lower.tail = x <= 0)) # if neg, lower tail, if pos, upper tail
    } else {
      return (pnorm(x, lower.tail = TRUE))
    }
  }
}

r2 <- function (y, x) sum(y * x)^2 / (sum(y^2) * sum(x^2))

## calculate partial correlation (the legit way) from dataframe ----

# Note now that this takes STRINGS so that covar_cols can take formula syntax
calc_pcor <- function(in_data, y_col, x_col, covar_cols) {
  stopifnot(is.character(y_col), is.character(x_col), is.character(covar_cols))
  if (length(covar_cols) > 1) {
    covar_cols <- paste(covar_cols, collapse = "+")
  }
  
  out <- in_data %>% 
    mutate(resid_y_covar = lm(as.formula(paste(y_col, covar_cols, sep = "~"))) %>% pluck("residuals"),
           resid_x_covar = lm(as.formula(paste(x_col, covar_cols, sep = "~"))) %>% pluck("residuals")) %>% 
    summarize(pcor = cor(resid_y_covar, resid_x_covar, method = "pearson")) %>% 
    pull(pcor)
  
  return (out)
}
