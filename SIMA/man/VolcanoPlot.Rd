% Generated by roxygen2 (4.0.1): do not edit by hand
\name{VolcanoPlot}
\alias{VolcanoPlot}
\title{Volcano plot}
\usage{
VolcanoPlot(fit, coef = NULL, FC = 2, PV = 0.05, p.val = c("adj",
  "raw"), name = "")
}
\arguments{
\item{fit}{an an ‘MArrayLM’ fitted linear model object.}

\item{coef,}{an Integer value corresponding to the desired contrast. Default to NULL.}

\item{FC,}{numeric, Fold Change threshold.}

\item{PV,}{numeric, p-value threshold.}

\item{p.val,}{character string, compute raw p-value ("raw") or adjusted p-value (adj"). Default to "adj". The p-values are adjusted with the procedure defined by Benjamini and Hochberg.}

\item{name,}{character string, title for the graphic. Default is an empty title ("").}
}
\description{
Volcano plot performed after a differential analysis. It represents the p-value obtained after a two sample (limma) test as a function of the Fold Change.
}

