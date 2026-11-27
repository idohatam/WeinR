

process_file_filter <- function(filepath, 
                                qc_obj, 
                                run_filter, 
                                min_avg_q, 
                                min_length,
                                max_n, 
                                adapter_sequence, 
                                max_mismatch_end, 
                                min_overlap_end,
                                min_internal_distance,
                                min_fragment_length,
                                start,
                                end,
                                base_dir,
                                output_filetype) {
  
  CHUNK_SIZE <- 10000L
  
  run_adapter <- !is.null(adapter_sequence)
  run_trim <- !is.null(start) || !is.null(end)
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  filename <- basename(filepath)
  single_qc <- qc_obj
  single_qc@files <- filepath
  single_qc@metrics[[filename]] <- qc_obj@metrics[[filename]]
  single_qc@summary_metrics <- qc_obj@summary_metrics[
    qc_obj@summary_metrics$File == filename, , drop = FALSE
  ]
  filter_summary <- new.env(parent = emptyenv())
  filter_summary$reads_before <- 0
  filter_summary$reads_after <- 0
  trim_summary <- new.env(parent = emptyenv())
  trim_summary$reads_before <- 0
  trim_summary$reads_after <- 0
  adapter_summary <- new.env(parent = emptyenv())
  adapter_summary$reads_before <- 0
  adapter_summary$reads_after <- 0
  
  message(paste("Processing:", filename))
  
  if (grepl("\\.bam$", filepath, ignore.case = TRUE)) {
    bam_file <- filepath
    bam_stream <- open(Rsamtools::BamFile(bam_file, yieldSize = CHUNK_SIZE))
    isBAM <- TRUE
    temp_sam <- paste0(base_dir, "/Processed_Files/", filename, ".sam")
    writeLines(c(
      "@HD\tVN:1.6\tSO:unsorted",
      "@SQ\tSN:*\tLN:0"
    ), con = temp_sam)
  } else {
    fastq_file <- filepath
    fastq_stream <- ShortRead::FastqStreamer(fastq_file, n = CHUNK_SIZE)
    isBAM <- FALSE
  }
  
  start_job <- Sys.time()
  start_time <- start_job
  
  start_index <- 1
  end_index <- CHUNK_SIZE
  
  repeat {
    
    if(isBAM) {
      chunk <- Rsamtools::scanBam(bam_stream)
      if(length(chunk[[1]][[1]]) == 0) {
        break
      }
      chunk_qsds <- QualityScaledDNAStringSet(chunk[[1]]$seq, chunk[[1]]$qual)
    } else {
      chunk <- ShortRead::yield(fastq_stream)
      if(length(chunk) == 0) {
        break
      }
      chunk_qsds <- as(chunk, "QualityScaledDNAStringSet")
    }
    if(length(chunk_qsds) != CHUNK_SIZE) {
      end_index <- start_index + length(chunk_qsds) - 1
    }
    # PROCESS CHUNK 
    if(run_filter) {
      processed_chunk <- filter_long(chunk_qsds, single_qc, filter_summary, start_index, end_index, min_avg_q, min_length, max_n)
    }
    if(run_adapter) {
      processed_chunk <- remove_adapter(processed_chunk, single_qc, adapter_summary, adapter_sequence)
    }
    if(run_trim) {
      processed_chunk <- trim_long(processed_chunk, single_qc, trim_summary, start_index, end_index, start, end)
    }
    start_index <- start_index + CHUNK_SIZE
    end_index <- end_index + CHUNK_SIZE
    output_file <- write_filter_output(processed_chunk, base_dir, filename, timestamp, output_filetype)
    # print(trim_summary$reads_before)
    # print(trim_summary$reads_after)
    
  }
  
  if(output_filetype == "bam") {
    Rsamtools::asBam(
      paste0(output_file, ".sam"),
      destination = output_file,
      overwrite   = TRUE
    )
  }
  
  
  
  if(isBAM) {
    close(bam_stream)
  } else {
    close(fastq_stream)
  }
  
  qc_obj@metadata[[filepath]]$filter_summary <- as.list(filter_summary)
  qc_obj@metadata[[filepath]]$adapter_summary <- as.list(adapter_summary)
  qc_obj@metadata[[filepath]]$trim_summary <- as.list(trim_summary)
  
  end_job <- Sys.time()
  print(end_job - start_job)
  message(paste0("Final output: ", output_file, ".", output_filetype))
  
  qc_obj
}