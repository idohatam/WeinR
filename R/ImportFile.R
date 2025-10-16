## A function that reads a sequence file
#  Args:
#    filePaths: a vector of strings, paths to the all the files
#  Returns:
#    -A list of QualityScaledDNAStringSet objects, a list item per file
#  Stops with an error if the file type isn't supported or file doesn't exist

ImportFile <- function(filePaths) {
  results <- lapply(filePaths, function(f) {
    infile <- fileType(f)
    if (infile %in% c("fastq", "fastq.gz")) {
      Biostrings::readQualityScaledDNAStringSet(f)
    } else {
      bam <- Rsamtools::scanBam(f)[[1]]
      Biostrings::QualityScaledDNAStringSet(bam$seq, bam$qual)
    }
  })
  names(results) <- basename(filePaths)
  return(results)
}