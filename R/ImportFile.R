## A function that reads a sequence file
#  Args:
#    filePaths: a vector of strings, paths to the all the files
#  Returns:
#    -A list of QualityScaledDNAStringSet objects, a list item per file
#  Stops with an error if the file type isn't supported or file doesn't exist

ImportFile <- function(filePath) {
    infile <- CheckFile(filePath)
    if (infile %in% c("fastq", "fastq.gz")) {
      Output <- Biostrings::readQualityScaledDNAStringSet(filePath)
    } else {
      bam <- Rsamtools::scanBam(filePath)[[1]]
      Output <- Biostrings::QualityScaledDNAStringSet(bam$seq, bam$qual)
    }
  
  result <- list(Output)
  names(result) <- filePath
  return(result)
}