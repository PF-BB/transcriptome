\name{clusterWithParameter}
\alias{clusterWithParameter}
\title{Draw a Cluster with Annotations}
\usage{
clusterWithParameter(X, pheno, title = " ", sample.name = colnames(X),
  dist.method = c("pearson", "euclidean", "manhattan", "spearman"),
  hclust.method = "ward", rect.hclust.k = 0)
}
\arguments{
  \item{X}{numeric matrix of the values to be plotted}

  \item{pheno}{a vector or matrix to annotate the samples}

  \item{title}{character string for the title of the
  cluster}

  \item{sample.name}{vector for names of the sample}

  \item{dist.method}{function used to compute the distance
  (dissimilarity) between both rows and columns}

  \item{hclust.method}{function used to compute the
  hierarchical clustering}

  \item{rect.hclust.k}{integer for numbers of groups : draw
  k rectangles}
}
\value{
\item{hclust object}{An object of class hclust.}
}
\description{
Cluster with annotations.
}

