\name{sima.heatmap}
\alias{sima.heatmap}
\title{Draw a Heat Map with Annotations}
\usage{
sima.heatmap(X, dist.method = c("pearson", "euclidean", "manhattan", "L1"),
  hclust.method = c("ward", "complete"), cols = c("RedGreen", "OrangeBlue",
  "RedBlue"), annot, labRow = rownames(X), labCol = colnames(X),
  title = "Heatmap", scale = "row", Rowv = TRUE, Colv = TRUE,
  hcr = NULL, hcc = NULL, saturate.cols = F)
}
\arguments{
  \item{X}{numeric matrix of the values to be plotted or
  object of class eSet}

  \item{dist.method}{function used to compute the distance
  (dissimilarity) between both rows and columns}

  \item{hclust.method}{function used to compute the
  hierarchical clustering}

  \item{cols}{colors of the heatmap}

  \item{annot}{a vector or matrix to annotate the samples}

  \item{labRow}{character vectors with row/column labels to
  use}

  \item{labCol}{character vectors with row/column labels to
  use}

  \item{title}{character string for the title of the
  cluster}

  \item{scale}{character indicating if the values should be
  centered and scaled ("row" or "none")}

  \item{Rowv}{booleen. TRUE plots samples/genes dendrogram.
  FALSE otherwise}

  \item{Colv}{booleen. TRUE plots samples/genes dendrogram.
  FALSE otherwise}

  \item{hcr}{object of class hclust if dendrogram already
  calculated}

  \item{hcc}{object of class hclust if dendrogram already
  calculated}

  \item{saturate.cols}{boolen.}
}
\value{
\item{list}{An list with 2 objects of class hclust.}
}
\description{
Heatmap with annotations.
}

