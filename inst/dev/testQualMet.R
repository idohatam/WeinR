folder_path <- "tests/testData"

# Get all .fastq files (non-recursive)
fastq_files <- list.files(path = folder_path, pattern = "\\.fastq$", full.names = TRUE)

# Print results
fastq_files

testObject <- .init_qc_object(fastq_files)

for(file in fastq_files){
  size_bytes <- file.info(file)$size
  size_mb <- size_bytes / 1024^2
  cat(
    "File:", file, "\n",
    "Size:", round(size_mb, 2), "MB\n"
  )
  start <- Sys.time()
  qs_list <- ImportFile(file)
 
  end <- Sys.time()
  print('file import')
  print(end - start)
  
  start <- Sys.time()
  testObject <-QualMat(testObject, qs_list, folder_path)
  end <- Sys.time()
  print('metrics calculation')
  print(end - start)
}


for(file in fastq_files){
  size_bytes <- file.info(file)$size
  size_mb <- size_bytes / 1024^2
  cat(
    "File:", file, "\n",
    "Size:", round(size_mb, 2), "MB\n"
  )
  start <- Sys.time()
  qs_list <- ImportFile(file)
  print('--------------')
  
  end <- Sys.time()
  print('file import')
  print(end - start)
  print('--------------')
  
  print('parallel per position calculation')
  start <- Sys.time()
  testchunked_quality_per_position(Biostrings::quality(qs_list))
  end <- Sys.time()
  print('overall')
  print(end - start)
  print('--------------')
  
  print('\nparallel + dynamic chunking per position calculation')
  start <- Sys.time()
  test2chunked_quality_per_position(Biostrings::quality(qs_list))
  end <- Sys.time()
  print('overall')
  print(end - start)
  print('--------------')
  
  start <- Sys.time()
  chunked_quality_per_position(Biostrings::quality(qs_list))
  end <- Sys.time()
  print('per position calculation')
  print(end - start)
  print('--------------')
  print('--------------')
  print('--------------')
  
}
