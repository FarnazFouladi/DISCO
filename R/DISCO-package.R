#' DISruption in COrrelations (DISCO)
#'
#' DISCO quantifies feature-specific changes in correlation patterns between
#' two groups, including changes in relationships across data modalities.
#'
#' @keywords internal
#' @importFrom dplyr %>% arrange bind_rows filter mutate
#' @importFrom stats na.omit p.adjust pcauchy pnorm
#' @importFrom tibble column_to_rownames
#' @importFrom utils write.csv
"_PACKAGE"

utils::globalVariables(c("abs_diff_cor", "p", "q", "r1", "r2"))
