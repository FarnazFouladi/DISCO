# Read historical cached analyses using the current presentation schema.
# Renaming and sorting only: no statistics are recomputed.
format_saved_disco_result <- function(result) {
  x <- result$disco
  old <- c(species = "Feature", num_corr = "Pairs", t = "DISCO_score",
           p = "DSICO_p", q = "DISCO_q")
  for (name in intersect(names(old), names(x))) names(x)[names(x) == name] <- old[[name]]
  columns <- c("Feature", "Pairs", "DISCO_score", "DSICO_p", "DISCO_q")
  x <- x[order(-x$DISCO_score, x$Feature), c(columns, setdiff(names(x), columns)), drop = FALSE]
  rownames(x) <- NULL
  result$disco <- x
  if (!is.null(result$cor_diff)) {
    names(result$cor_diff)[names(result$cor_diff) == "t1"] <- "feature_1"
    names(result$cor_diff)[names(result$cor_diff) == "t2"] <- "feature_2"
  }
  if (!is.null(result$feature_null_summary))
    names(result$feature_null_summary)[names(result$feature_null_summary) == "species"] <- "Feature"
  result
}
read_disco_result <- function(path) format_saved_disco_result(readRDS(path))
