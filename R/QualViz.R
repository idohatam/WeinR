#'QualViz is a one stop shop facilitator function 
#'


QualViz <- function(files){
  
  qc_obj <-.init_qc_object(files)
  
  for(file in qc_obj@files){
    
    #Import file, return as a list of QualityScaledDNAStringSet
    qstringset <- ImportFile(file)
    
    #Calculate all the quality stuff
    qc_obj <- QualMat(qc_obj, qstringset, file)
    
    #Plot quality plots
    #qc_obj <- QualPlot(qc_obj)
  }
  
  #Report generation
  #ReportFun(qc_obj)
  
  #Return the LongReadQC object
  base::return(qc_obj)
  
}
