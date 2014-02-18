\name{acplight}
\alias{acplight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
acplight(titre = "", logValues, classes, couleurs, formes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{titre}{
%%     ~~Describe \code{titre} here~~
}
  \item{logValues}{
%%     ~~Describe \code{logValues} here~~
}
  \item{classes}{
%%     ~~Describe \code{classes} here~~
}
  \item{couleurs}{
%%     ~~Describe \code{couleurs} here~~
}
  \item{formes}{
%%     ~~Describe \code{formes} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (titre = "", logValues, classes, couleurs, formes = NULL) 
{
    if (is.null(formes)) {
        formes <- 16
        ff <- 16
    }
    else {
        ff <- as.numeric(levels(as.factor(formes)))
    }
    PRC <- prcomp(t(logValues), scale = TRUE)
    PC_prc <- PRC$x
    cc <- palette()[as.numeric(levels(as.factor(couleurs)))]
    prcplot <- xyplot(PC_prc[, 2] ~ PC_prc[, 1], groups = couleurs, 
        col = cc, xlab = "PC1", ylab = "PC2", main = titre, pch = ff, 
        key = list(space = "right", points = list(tab = classes, 
            pch = ff), text = list(tab = classes, pch = ff), 
            col = cc))
    plot(prcplot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
