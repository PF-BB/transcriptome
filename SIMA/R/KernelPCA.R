KernelPCA <-
function(X) {
	K   <- X%*%t(X)
	l <- dim(X)[1]
	# Centrer
	D <- apply(K,2,mean)
	E <- mean(D)
	j <- matrix(1,l,1)%*%D
	K <- K-j-t(j)+E*matrix(1,l,l)
	
	# # Reduire
	# D <- diag(1/(sqrt(diag(K))))
	# K <- D %*% K %*% D
	
	vec <- eigen(K,symmetric=T)$vectors
	vp  <- abs(eigen(K,symmetric=T)$values)
	vpmat <- t(matrix(rep(vp,l),l,l))
	alpha <- 1/sqrt(vpmat)*vec
	PC <- K%*%alpha

	list(PC=PC)
}
