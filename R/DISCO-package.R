#' DISruption in COrrelations (DISCO)
#'
#' DISCO quantifies feature-specific changes in correlation patterns between
#' two groups, including changes in relationships across data modalities.
#'
#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom stats pnorm
#' @importFrom tibble column_to_rownames
"_PACKAGE"

utils::globalVariables(c(
  "filter_threshold",
  "annotation",
  "lower",
  "percentage",
  "neg_log10_real_disco_p",
  "null_discovery_percentage",
  "panel",
  "real_q_category",
  "real_to_null_ratio",
  "result_type",
  "upper"
))
