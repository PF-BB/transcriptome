\name{VolcanoPlot}
\alias{VolcanoPlot}
\title{Volcano plot}
\usage{
VolcanoPlot(topT, FC, PV, p.val = c("adj", "raw"), name = "", feature.id)
}
\arguments{
  \item{topT}{object topTable}

  \item{FC}{numeric. Seuil de Fold change}

  \item{PV}{numeric. Seuil de p-value}

  \item{p.val}{compute of raw pvalue or on adjusted pvalue}

  \item{name}{title for the graphic}
}
\description{
Volcano plot.
}

