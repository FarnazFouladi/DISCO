# Write tabular outputs without an extra row-number column.
write_result_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}
