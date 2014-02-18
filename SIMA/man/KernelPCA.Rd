\name{KernelPCA}
\alias{KernelPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KernelPCA(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
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
function (X) 
{
    K <- X \%*\% t(X)
    l <- dim(X)[1]
    D <- apply(K, 2, mean)
    E <- mean(D)
    j <- matrix(1, l, 1) \%*\% D
    K <- K - j - t(j) + E * matrix(1, l, l)
    vec <- eigen(K, symmetric = T)$vectors
    vp <- abs(eigen(K, symmetric = T)$values)
    vpmat <- t(matrix(rep(vp, l), l, l))
    alpha <- 1/sqrt(vpmat) * vec
    PC <- K \%*\% alpha
    list(PC = PC)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
