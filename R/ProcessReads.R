

process_reads <- function(qc_obj,
                          filter = FALSE,
                          MinAvgQS = 5,
                          MinLength = 100,
                          MaxNumberNs = 2,
                          AdapterSeq = NULL,
                          MaxMismatchEnd = 3,
                          MinOverlapEnd = 20,
                          MinInternalDistance = 100,
                          MinFragmentLength = 200,
                          Start = NULL,
                          End = NULL,
                          OutFileType = c("fastq"),
                          base_dir = "./WeinR_Outputs") {
  
  message("Start processing ", length(qc_obj@files), " file(s)...")
  
  for (file in qc_obj@files) {
    qc_obj <- process_file_filter(file, qc_obj, filter, MinAvgQS, MinLength, MaxNumberNs, AdapterSeq, MaxMismatchEnd,
                                  MinOverlapEnd, MinInternalDistance, MinFragmentLength, Start, End, base_dir, OutFileType)
  }
  
  qc_obj
}