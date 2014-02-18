\name{scalefactor}
\alias{scalefactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scalefactor(noms, tv = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{noms}{
%%     ~~Describe \code{noms} here~~
}
  \item{tv}{
%%     ~~Describe \code{tv} here~~
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
function (noms, tv = 500) 
{
    mat <- sapply(noms, function(x) exprs(mas5(ReadAffy(filenames = x), 
        normalize = FALSE)))
    colnames(mat) <- noms
    sf <- apply(mat, 2, function(x) tv/mean(x, trim = 0.02))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
