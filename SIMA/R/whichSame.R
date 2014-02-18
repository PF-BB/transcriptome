whichSame <-
function(listeVrac,listeRef) {
	foo     <- function(x,listChain)  grep( x , listChain )
	indices <- unlist(sapply(listeVrac,foo,listeRef))
}
