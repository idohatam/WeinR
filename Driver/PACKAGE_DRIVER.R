
# PACKAGE DRIVER
#change this as well
working_directory <- "/Users/aaronsantos/BIO4415/WeinR"
setwd(working_directory)

function_names <- list.files("./R")
setwd("./R")
lapply(function_names, source)
setwd("../")

GENUS <- "Mus"

#change file path 

qc_obj <- create_report(files = "/Users/aaronsantos/Downloads/subsampled_1percent.bam", report_name = "GENUS", render_report = TRUE, force = TRUE)



