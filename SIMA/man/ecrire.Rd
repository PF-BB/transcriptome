\name{ecrire}
\alias{ecrire}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ecrire(logValues, top, titre = "Tableau.tsv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logValues}{
%%     ~~Describe \code{logValues} here~~
}
  \item{top}{
%%     ~~Describe \code{top} here~~
}
  \item{titre}{
%%     ~~Describe \code{titre} here~~
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
function (logValues, top, titre = "Tableau.tsv") 
{
    Ncol <- dim(logValues)[2]
    g <- dim(logValues)[1]
    Ncomp <- length(top)
    print(c(Ncol, g, Ncomp))
    tableau <- as.data.frame(matrix(0, g, 1 + Ncol + 4 * Ncomp))
    ordreID <- rownames(logValues)
    noms_tab <- vector(length = 1 + Ncol + 4 * Ncomp)
    tableau[, 1] <- ordreID
    noms_tab[1] <- "ID"
    tableau[, 1 + 1:Ncol] <- logValues
    noms_tab[1 + 1:Ncol] <- colnames(logGCRMA)
    for (ind in 1:Ncomp) {
        print(dim(tableau))
        temp <- top[[ind]]
        print(names(temp))
        temp <- temp[order(temp$ID), ]
        DF <- sign(temp$logFC) * 2^abs(temp$logFC)
        tableau[, 1 + Ncol + 4 * (ind - 1) + 1:4] <- cbind(temp$logFC, 
            DF, temp$P.Value, temp$adj.P.Val)
        noms_tab[1 + Ncol + 4 * (ind - 1) + 1:4] <- paste(names(top[ind]), 
            c("logFC", "DF", "p-value", "adjusted p-value"))
    }
    colnames(tableau) <- noms_tab
    print(tableau[1:10, ])
    write.table(tableau, titre, dec = ",", sep = "\t", row.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
