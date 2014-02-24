#' Declare the contrasts of interests for the differential analysis.
#' @param y A vector containing the class of each sample.
#' @return \item{contr}{The vector containing the character strings describing the desired constrasts.}
#' @title Declare contrasts.
#' @export declare_contrastes

declare_contrastes <- function(y){
  
  # 1. Check y
  if (missing(y)) stop("\n\tL'argument y est manquant !")
  if (!is.factor(y)) warning("y n'est pas un facteur, il sera transformé automatiquement en facteur.")

  yfac <- as.factor(y)
  lev <- levels(yfac)
  nlev <- length(lev)
  
  # 2. Create the contrasts
  
  combi <- t( combn(lev, 2) )
  ncomb <- nrow(combi)
  MC <- matrix(0, nlev, ncomb, dimnames =  list(lev, paste0("Cont",1:ncomb) )  )
  for (i in 1:ncomb) MC[combi[i,],i] <- c(1,-1)
  
  MC <- fix(MC)
  contr <- apply( MC, 2, function(v) if(all(v==0)) {
                                        return(character(0)) 
                                      } else if ( !identical(names(table(v)), c("-1","0","1")) ) {
                                        stop("Utiliser seulement -1, 0 or 1 pour declarer des contrastes.")
                                      } else if ( (sum(v==1) != 1) | (sum(v==-1)!= 1) ) {
                                        stop("Trop de 1 ou de -1 pour un contraste.")
                                      } else {
                                        return( sprintf("%s - %s", lev[v == 1], lev[v == -1])) 
                                      } )  
  unlist(contr)
  return(unlist(contr))
  
}



