#'QualViz is a one stop shop facilitator function 
#'


QualViz <- function(files){
  
  qc_obj <-.init_qc_object(files)
  
  #Check files returns a vector with either bam, fastq, or mixed
  file_type <- FileCheck(files = qc_obj@files)
  
  #Import files, returnes them as a list of QualityScaledDNAStringSet
  
  qstringsets <- ImoportFiles(files = qc_obj@files, ftype = file_type)
  
  #Calc all the quality stuff
  
  qc_obj <- QualMat(qc_obj)
  
  #Plot quality plots
  
  qc_obj <- QualPlot(qc_obj)
  
  #Report generation
  
  ReportFun(qc_obj)
  
  #Return the LongReadQC object
  base::return(qc_obj)
  
}