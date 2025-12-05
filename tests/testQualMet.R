folder_path <- "/Users/navameadi/Downloads/Nakazawaea_holstii.fastq"

# Get all .fastq files (non-recursive)
fastq_files <- list.files(path = folder_path, pattern = "\\p.fastq$", full.names = TRUE)

# Print results
fastq_files

testObject <- .init_qc_object(folder_path)

for(file in folder_path){
  # start <- Sys.time()
  qs_list <- ImportFile(file)
  testObject <-QualMat(testObject, qs_list, folder_path)
  # end <- Sys.time()
  # print('file import')
  # print(end - start)
  
  # start <- Sys.time()
  # testObject <- BenchmarkQualMat(testObject, qs_list, file)
  # # end <- Sys.time()
  # print('metrics calculation')
  # print(end - start)
}
