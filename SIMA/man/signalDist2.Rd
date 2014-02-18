\name{signalDist2}
\alias{signalDist2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
signalDist2(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
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
function (object) 
{
    par(mfrow = c(2, 1))
    ArrayIndex = as.character(1:length(sampleNames(object)))
    boxplot(object, names = ArrayIndex, ylab = "Log2(Intensity)", 
        xlab = "Array Index")
    hist(x = object, lt = 1:length(ArrayIndex), col = 1:length(ArrayIndex), 
        xlab = "Log2(Intensity)", which = "both")
    temppar <- par()
    legend(((temppar$xaxp[2] - temppar$xaxp[1])/temppar$xaxp[3]) * 
        (temppar$xaxp[3] - 1) + temppar$xaxp[1], temppar$yaxp[2], 
        as.character(ArrayIndex), lt = 1:length(ArrayIndex), 
        col = 1:length(ArrayIndex), cex = 0.5)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
