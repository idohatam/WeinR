## A function that reads a sequence file
#  Args:
#    filePaths: a vector of strings, paths to the all the files
#  Returns:
#    -A list of QualityScaledDNAStringSet objects, a list item per file
#  Stops with an error if the file type isn't supported or file doesn't exist

## A function that reads a sequence file
#  Args:
#    filePaths: a vector of strings, paths to the all the files
#  Returns:
#    -A list of QualityScaledDNAStringSet objects, a list item per file
#  Stops with an error if the file type isn't supported or file doesn't exist

ImportFile <- function(filePath) {
  infile <- CheckFile(filePath)
  
  if (infile %in% c("fastq", "fastq.gz")) {
    # For FASTQ files, Biostrings handles names automatically
    Output <- Biostrings::readQualityScaledDNAStringSet(filePath)
    
  } else if (infile == "bam") {
    # For BAM files, use scanBam to get sequences, qualities, and names
    bam <- Rsamtools::scanBam(filePath)[[1]]
    
    # Create a QualityScaledDNAStringSet
    Output <- Biostrings::QualityScaledDNAStringSet(bam$seq, bam$qual)
    
    # Attach read names (qname) if available
    if (!is.null(bam$qname) && length(bam$qname) == length(Output)) {
      names(Output) <- bam$qname
    }
  } 
  
  result <- Output
  return(result)
}
