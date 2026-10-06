# Bound p-values before applying the existing Cauchy transformation
safe_tan_cauchy <- function(p_vect, tol = 1e-16) {
  p_vect <- pmax(pmin(p_vect, 1 - tol), tol)
  1 / length(p_vect) * sum(tan((0.5 - p_vect) * pi))
}
