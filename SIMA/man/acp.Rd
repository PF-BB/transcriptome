\name{acp}
\alias{acp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
acp(fichier, logValues, classes, couleurs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fichier}{
%%     ~~Describe \code{fichier} here~~
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
function (fichier, logValues, classes, couleurs) 
{
    PRC <- prcomp(t(logValues), scale = TRUE)
    PC_prc <- PRC$x
    noms <- colnames(logValues)
    cc <- palette()[as.numeric(levels(as.factor(couleurs)))]
    cc2 <- palette()[couleurs]
    pourcent <- xyplot(100 * (summary(PRC)$importance[3, ]) ~ 
        1:length(PRC$sdev), main = "Cumulative proportion of Variance", 
        xlab = "PC number", ylab = "Percentage", panel = function(x, 
            y) {
            llines(x, y)
            lpoints(x, y, pch = 5)
        })
    prcplot <- xyplot(PC_prc[, 2] ~ PC_prc[, 1], groups = couleurs, 
        col = cc, xlab = "PC1", ylab = "PC2", main = "2D-PCA", 
        pch = 16, key = list(space = "right", points = list(tab = classes, 
            pch = 16), text = list(tab = classes, pch = 16), 
            col = cc))
    prctextplot <- xyplot(PC_prc[, 2] ~ PC_prc[, 1], groups = noms, 
        panel = function(x, y, subscripts, groups) ltext(x = x, 
            y = y, label = groups[subscripts], cex = 1, fontfamily = "HersheySans", 
            col = cc2))
    pdf(fichier, sep = "", collapse = "")
    plot(prcplot)
    biplot(PRC)
    plot(update(prctextplot, aspect = "xy", xlab = "PC1", ylab = "PC2", 
        main = "2D-PCA"))
    plot(pourcent)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
