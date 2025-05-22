cohens.d <- function (vec1, vec0) {
  (mean(vec1) - mean(vec0)) / sqrt((var(vec1)+var(vec0))/2)
}
