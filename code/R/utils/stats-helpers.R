cohens.d <- function (vec1, vec0, na.rm = FALSE) {
  (mean(vec1, na.rm = na.rm) - mean(vec0, na.rm = na.rm)) / sqrt((var(vec1, na.rm = na.rm)+var(vec0, na.rm = na.rm))/2)
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
