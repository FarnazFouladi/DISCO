# Apply the equal-weight Cauchy transformation using the nearest p-value tail
safe_tan_cauchy <- function(p_vect) {
  tail_boundary <- 1e-8
  lower <- p_vect <= 0.5
  tail_p <- pmax(pmin(p_vect, 1 - p_vect), tail_boundary)
  terms <- ifelse(lower, 1, -1) / tan(pi * tail_p)
  terms[p_vect == 0.5] <- 0

  mean(terms)
}
