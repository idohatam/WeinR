folder_path <- "/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/Sungouiella/"

# Get all .fastq files (non-recursive)
fastq_files <- list.files(path = folder_path, pattern = "\\p.fastq$", full.names = TRUE)

# Print results
fastq_files

testObject2 <- .init_qc_object(fastq_files)

for(file in fastq_files){
  start <- Sys.time()
  qs_list <- ImportFile(file)
  end <- Sys.time()
  print('file import')
  print(end - start)
  
  start <- Sys.time()
  testObject2 <- BenchmarkQualMat(testObject2, qs_list, file)
  end <- Sys.time()
  print('metrics calculation')
  print(end - start)
}