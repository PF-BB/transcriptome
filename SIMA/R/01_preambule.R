#' Create hierarchy
#' @param path where to create the folders
#' @references 
#' @export preambule


preambule <- function(path_in="./",path_out="./"){
  
  if(!file.exists(path_in))
    stop("Le dossier ", path_in," n'existe pas ! ", call.=FALSE)
  else if(!file.exists(path_out))
    stop("Le dossier ", path_out," n'existe pas ! ", call.=FALSE)
  else if(!file.exists(file.path(path_in,"phenotype.txt")))
    stop("Le fichiers de donnee phenotypiques ", file.path(path_in,"phenotype.txt")," n'existe pas ! ", call.=FALSE)
  else{
    cat("Création de l'arborescence ... \n\n")
    if(!file.exists(paste(path,"QC",sep="/")))
      dir.create(paste(path,"QC",sep="/"))
    
    if(!file.exists(paste(path,"AnalyseDifferentielle",sep="/")))
      dir.create(paste(path,"AnalyseDifferentielle",sep="/"))

    if(!file.exists(paste(path,"Data",sep="/")))
      dir.create(paste(path,"Data",sep="/"))
    
    if(!file.exists(paste(path,"Graphes",sep="/")))
      dir.create(paste(path,"Graphes",sep="/"))
    
  }
  
  
}