#' Create an output arborescence in which the output files will be saved. 
#' @param path_in the path of the folder containing the input data.
#' @param path_out the path of the folder where the ouput files should be written.
#' @export preambule


preambule <- function(path_in="./",path_out="./"){
  
  if(!file.exists(path_in))
    stop("Le dossier ", path_in," n'existe pas ! ", call.=FALSE)
  else if(!file.exists(path_out))
    stop("Le dossier ", path_out," n'existe pas ! ", call.=FALSE)
  else if(!file.exists(file.path(path_in,"phenotype.txt")))
    stop("Le fichiers de donnee phenotypiques ", file.path(path_in,"phenotype.txt")," n'existe pas ! ", call.=FALSE)
  else{
    cat("Creation de l'arborescence ... \n\n")
    if(!file.exists(paste(path_out,"QC",sep="/")))
      dir.create(paste(path_out,"QC",sep="/"))
    
    if(!file.exists(paste(path_out,"AnalyseDifferentielle",sep="/")))
      dir.create(paste(path_out,"AnalyseDifferentielle",sep="/"))

    if(!file.exists(paste(path_out,"Data",sep="/")))
      dir.create(paste(path_out,"Data",sep="/"))
    
    if(!file.exists(paste(path_out,"Graphes",sep="/")))
      dir.create(paste(path_out,"Graphes",sep="/"))
    
  }
  
  
}
