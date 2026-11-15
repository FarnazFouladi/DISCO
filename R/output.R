# Excel permits 1,048,576 rows per sheet, including the column headings.
# Split large pairwise tables across sheets rather than losing rows.
write_result_xlsx <- function(x, path, max_rows = 1048575L) {
  starts <- seq.int(1L, max(1L, nrow(x)), by = max_rows)
  sheets <- lapply(starts, function(start) {
    if (!nrow(x)) return(x)
    x[seq.int(start, min(start + max_rows - 1L, nrow(x))), , drop = FALSE]
  })
  names(sheets) <- if (length(sheets) == 1L) "Results" else paste0("Results_", seq_along(sheets))
  writexl::write_xlsx(sheets, path)
  invisible(path)
}
