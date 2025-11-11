folder_path <- "/Users/ocoelho00/Documents/biol4315/repos/longR/testData"

# Get all .fastq files (non-recursive)
fastq_files <- list.files(path = folder_path, pattern = "\\_1.fastq$", full.names = TRUE)

# Print results
fastq_files

qs_list <- ImportFile(fastq_files)

testObject <- .init_qc_object(fastq_files)

testObjectMetricd <- QualMat(testObject, qs_list)

