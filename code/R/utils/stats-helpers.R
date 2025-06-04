cohens.d <- function (vec1, vec0) {
  (mean(vec1) - mean(vec0)) / sqrt((var(vec1)+var(vec0))/2)
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
