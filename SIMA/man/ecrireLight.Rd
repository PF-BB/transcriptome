\name{ecrireLight}
\alias{ecrireLight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ecrireLight(logValues, top, fichier = "tableau.tsv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logValues}{
%%     ~~Describe \code{logValues} here~~
}
  \item{top}{
%%     ~~Describe \code{top} here~~
}
  \item{fichier}{
%%     ~~Describe \code{fichier} here~~
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
function (logValues, top, fichier = "tableau.tsv") 
{
    smbls <- as.list(hgu133plus2SYMBOL)
    Ncol <- dim(logValues)[2]
    g <- dim(logValues)[1]
    ordreID <- rownames(logValues)
    temp <- top
    print(names(temp))
    temp <- temp[order(temp$ID), ]
    DF <- sign(temp$logFC) * 2^abs(temp$logFC)
    tableau <- dataframe(ID = ordreID, logValues, logFC = temp$logFC, 
        DF = DF, `p-value` = temp$P.Value, `adjusted p-value` = temp$adj.P.Val)
    print(tableau[1:10, ])
    write.table(tableau, titre, dec = ",", sep = "\t", row.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
