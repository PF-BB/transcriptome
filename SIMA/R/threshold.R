threshold <-
function(vec,seuil){
  out <- rep(0,length(vec))
  out[which(vec > seuil)] <- 1
  out
}
