\name{QCReport_simple}
\alias{QCReport_simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QCReport_simple(projet, nom_Etude)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projet}{
%%     ~~Describe \code{projet} here~~
}
  \item{nom_Etude}{
%%     ~~Describe \code{nom_Etude} here~~
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
function (projet, nom_Etude) 
{
    QCproj <- qc(projet)
    ind1 <- 1
    ind2 <- 3
    ind3 <- 4
    ind4 <- 6
    prob1 <- substr(colnames(QCproj@qc.probes)[ind1], 1, nchar(colnames(QCproj@qc.probes)[ind1]) - 
        5)
    prob2 <- substr(colnames(QCproj@qc.probes)[ind3], 1, nchar(colnames(QCproj@qc.probes)[ind3]) - 
        5)
    ratios <- cbind(QCproj@qc.probes[, ind1] - QCproj@qc.probes[, 
        ind2], QCproj@qc.probes[, ind3] - QCproj@qc.probes[, 
        ind4])
    rownames(ratios) <- rownames(QCproj@qc.probes)
    colnames(ratios) <- c(prob1, prob2)
    expratios <- 2^ratios
    pdf(paste(c(nom_Etude, "_QC.pdf"), sep = "", collapse = ""))
    afficheTable3col(2^ratios)
    plotQC2(QCproj, cex = 0.5)
    plot.new()
    titlePage(projet)
    signalDist2(projet)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
