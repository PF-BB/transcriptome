\name{topTableReducer}
\alias{topTableReducer}
\title{Top table reducer.}
\usage{
topTableReducer(fit, coef = NULL, FC = 2, PV = 0.05, p.val = c("adj",
  "raw"), tlog2 = FALSE)
}
\arguments{
  \item{fit}{.}

  \item{FC}{.}

  \item{PV}{.}

  \item{p.val}{.}
}
\value{
\item{top}{An object topTable with a twist obtained after
using limma's fit function.}
}
\description{
Modify the original topTable function to allow non-log fold
changes and average expressions.
}
\references{
limma.
}

