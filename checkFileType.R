## A function that checks if a file exists and determines its type
#  Args:
#    filePath: string, path to the file
#  Returns:
#    Type of file ("bam", "fastq", or "fastq.gz")
#  Stops with an error if the file doesn't exist or the file type isn't compatible

fileType <- function(filePath) 
  {
  # check if file exists, if not stop with an error message
  if(!file.exists(filePath))
  {
    stop(paste("File does not exist:", filePath))
  }
  if (grepl("\\.bam$", filePath, ignore.case = TRUE)) 
  {
    return("bam")
  } 
  else if (grepl("\\.(fastq|fq)\\.gz$", filePath, ignore.case = TRUE)) 
  {
    return("fastq.gz")
  } 
  else if (grepl("\\.(fastq|fq)$", filePath, ignore.case = TRUE)) 
  {
    return("fastq")
  } 
  else 
  {
    # if file isn't supported, stop with an error message
    stop("Usupported file type. Expected bam or fastq file.")
  }
}
