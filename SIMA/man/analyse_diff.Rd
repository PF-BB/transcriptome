\name{analyse_diff}
\alias{analyse_diff}
\title{Perform differential analysis.}
\usage{
analyse_diff(X, y, contr)
}
\arguments{
  \item{X}{An Expression Set.}

  \item{y}{A vector contaning the class whose each sample
  blongs to.}

  \item{contr}{A vector of character strongs describing the
  contrasts.}
}
\value{
\item{top}{An object topTable obtained after using limma's
fit function.}
}
\description{
Declare the contrasts of interests for the differential
analysis.
}
\references{
limma.
}

