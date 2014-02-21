#' Heatmap with annotations.
#' @param eSet an expression set
#' @references 
#' @title Draw a density plot for an eSet
#' @export sima.densityplot


sima.densityplot <- function (eSet) {
  plot(eSet, what="density")
}
