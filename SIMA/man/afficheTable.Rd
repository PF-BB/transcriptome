\name{afficheTable}
\alias{afficheTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
afficheTable(object)
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
    polygon(c(0, 0, 0.9, 0.9, 0), c(0.05, 0.95, 0.95, 0.05, 0.05))
    polygon(c(0.45, 0.45), c(0.05, 0.95))
    polygon(c(0, 0.9), c(0.88, 0.88))
    polygon(c(0, 0.9), c(0.87, 0.87))
    text(0.5, 1, "Table des ratios 3' sur 5'", cex = 1.5)
    text(0.4, 0.03, date(), cex = 0.5)
    text(0.45/2, 0.9, colnames(object)[1], cex = 0.8)
    text(3 * 0.45/2, 0.9, colnames(object)[2], cex = 0.8)
    n <- dim(object)[1]
    cexval <- min(10/n, 1)
    for (i in 1:n) {
        text(0.45/2, 0.9 - i * (0.9 - 0.01)/(n + 1), round((object)[i, 
            1], 2), cex = cexval)
        polygon(c(0, 0.9), c(0.9 - (i + 0.5) * (0.9 - 0.01)/(n + 
            1), 0.9 - (i + 0.5) * (0.9 - 0.01)/(n + 1)))
        text(3 * 0.45/2, 0.9 - i * (0.9 - 0.01)/(n + 1), round((object)[i, 
            2], 2), cex = cexval)
    }
    return(TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
