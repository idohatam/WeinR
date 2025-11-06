#'QualViz is a one stop shop facilitator function 
#'


QualViz <- function(files){
  
  qc_obj <-.init_qc_object(files)
  
  for(i in seq_along(qc_obj@files)){
    
    #Import file, return as a list of QualityScaledDNAStringSet
    qstringset <- ImoportFiles(qc_obj@file$i)
    
    #Calculate all the quality stuff
    qc_obj <- QualMat(qc_obj, qstringset, qc_obj@files$i)
    
    #Plot quality plots
    #qc_obj <- QualPlot(qc_obj)
  }
  
  #Report generation
  #ReportFun(qc_obj)
  
  #Return the LongReadQC object
  base::return(qc_obj)
  
}