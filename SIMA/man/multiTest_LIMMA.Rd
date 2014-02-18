\name{multiTest_LIMMA}
\alias{multiTest_LIMMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
multiTest_LIMMA(X, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
function (X, y) 
{
    g <- nrow(X)
    design <- cbind(y, 1 - y)
    colnames(design) <- c("un", "deux")
    cont <- makeContrasts(un - deux, levels = design)
    fit <- lmFit(X, design)
    fit <- contrasts.fit(fit, cont)
    ordinary.t <- fit$coef/fit$stdev.unscaled/fit$sigma
    fit <- eBayes(fit)
    top <- topTable(fit, number = g, adjust.method = "BH", sort.by = "P")
    list(top = top, ordinary.t = ordinary.t)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
