\name{lumi_normalization}
\alias{lumi_normalization}
\title{Normalization of Illumina}
\usage{
lumi_normalization(dataFile, pheno, bg.method = "none",
  norm.method = "quantile")
}
\arguments{
  \item{dataFile}{A character string for BeadStudio file.}

  \item{pheno}{An oject of class AnnotatedDataFrame.}

  \item{bg.method}{A character string: the method used to
  deals with microarray background noise. If the specified
  method is not supported, an error is generated}

  \item{norm.method}{A character string: the method used
  perform an inter-array normalization. If the specified
  method is not supported, an error is generated}
}
\value{
\item{eset}{An object of classe Expression-set.}
}
\description{
Normalization of Illumina datasets.
}
\references{
lumi.
}

