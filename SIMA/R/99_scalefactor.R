scalefactor <-
function (noms, tv = 500){
  mat <- sapply(noms, function(x) exprs(mas5(ReadAffy(filenames=x),normalize=FALSE)) )
  colnames(mat) <- noms
  sf <- apply(mat,2,function(x) tv/mean(x,trim=0.02))
}
