RemoveExt <- function(x) {
  sub("\\.(fastq|fq|fasta|fa|bam)(\\.gz)?$", "", basename(x), ignore.case = TRUE)
}