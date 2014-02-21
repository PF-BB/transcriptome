#' Create hierarchy
#' @param path where to create the folders
#' @references 
#' @export preambule


preambule <- function(path="./"){
  
  if(!file.exists(path))
    stop("Le dossier ", path," n'existe pas ! ", call.=FALSE)
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