write_fastq <- function(path, records) {
  stopifnot(is.character(path), length(path) == 1L)
  stopifnot(is.list(records), length(records) >= 1L)
  
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  for (rec in records) {
    stopifnot(all(c("id","seq","plus","qual") %in% names(rec)))
    writeLines(c(
      paste0("@", rec$id),
      rec$seq,
      rec$plus,
      rec$qual
    ), con = con)
  }
  invisible(path)
}

write_fastq_gz <- function(path, records) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  for (rec in records) {
    writeLines(c(
      paste0("@", rec$id),
      rec$seq,
      rec$plus,
      rec$qual
    ), con = con)
  }
  invisible(path)
}

minimal_records_ok <- function(n = 2L) {
  lapply(seq_len(n), function(i) {
    list(
      id = paste0("read", i),
      seq = "ACGTACGT",
      plus = "+",
      qual = "FFFFFFFF"
    )
  })
}


