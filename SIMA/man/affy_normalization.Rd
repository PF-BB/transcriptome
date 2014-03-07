\name{affy_normalization}
\alias{affy_normalization}
\title{Normalization of Affymetrix}
\usage{
affy_normalization(inputFolder, pheno, bg.method = "none",
  norm.method = "quantile")
}
\arguments{
  \item{inputFolder}{A character string: where to find the
  CEL files and the phenotype data.}

  \item{pheno}{An oject of class AnnotatedDataFrame.}

  \item{bg.method}{A character string: the method used to
  deals with microarray background noise. If the specified
  method is not supported, an error is generated}

  \item{norm.method}{A character string: the method used
  perform an inter-array normalization. If the specified
  method is not supported, an error is generated}
}
\value{
\item{eset}{An object of class Expression-set.}
}
\description{
Normalization of Affymetrix datasets.
}
\references{
FARMS.
}

