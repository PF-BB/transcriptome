\name{normalization}
\alias{normalization}
\title{Normalization of Affymetrix or Illumina datasets.}
\usage{
normalization(input, phenoFile, type, bg.method = "none",
  norm.method = "quantile", bplot = TRUE, bwrite = TRUE, output = "./")
}
\arguments{
  \item{input}{A character string: either the name of the
  folder where to find the CEL files and the phenotype data
  OR the name of the raw Illumina file obtained with
  BeadStudio.}

  \item{phenoFile}{A character string: where to find the
  phenotypes of interest. It is required to be a
  tab-separated file.}

  \item{type}{A character string: must be either "affy" or
  "illumina".}

  \item{bg.method}{A character string: the method used to
  deals with microarray background noise. If the specified
  method is not supported, an error is generated}

  \item{norm.method}{A character string: the method used
  perform an inter-array normalization. If the specified
  method is not supported, an error is generated}

  \item{write}{A boolean to indicate if the data should be
  written}
}
\value{
\item{eset}{An object of classe Expression-set.}
}
\description{
Normalization of Affymetrix or Illumina datasets.
}
\references{
FARMS and lumi.
}

