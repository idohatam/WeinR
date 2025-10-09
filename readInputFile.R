## A function that reads a sequence file
#  Args:
#    filePaths: a vector of strings, paths to the all the files
#  Returns:
#    -A list of QualityScaledDNAStringSet objects, a list item per file
#  Stops with an error if the file type isn't supported or file doesn't exist

readInputFiles <- function(filePaths) 
  {
  lapply(filePaths, function(f) 
    {
    infile <- fileType(f)
    basename <- basename(f)
    
    if (infile %in% c("fastq", "fastq.gz")) 
    {
      result <- Biostrings::readQualityScaledDNAStringSet(f)
      
    } 
    else
    {
       bam <-Rsamtools::scanBam(f)[[1]]
       result <- Biostrings::QualityScaledDNAStringSet(bam$seq, bam$qual)
    } 
    setNames(list(result), basename)
  })
}

