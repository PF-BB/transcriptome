\name{volcano}
\alias{volcano}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
volcano(top, seuilFC, seuilPval, correction = correction, titre = "Volcano Plot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{top}{
%%     ~~Describe \code{top} here~~
}
  \item{seuilFC}{
%%     ~~Describe \code{seuilFC} here~~
}
  \item{seuilPval}{
%%     ~~Describe \code{seuilPval} here~~
}
  \item{correction}{
%%     ~~Describe \code{correction} here~~
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
function (top, seuilFC, seuilPval, correction = correction, titre = "Volcano Plot") 
{
    if (correction) {
        B <- -log10(top$adj.P.Val)
        lab <- "-log10(adjusted p-value)"
        titre <- paste(c(titre, " (FDR)"), collapse = "", sep = "")
    }
    else {
        B <- -log10(top$P.Value)
        lab <- "-log10(p-value)"
    }
    FC <- top$logFC
    gauche <- intersect(which(B > -log10(seuilPval)), which(FC < 
        -seuilFC))
    droite <- intersect(which(B > -log10(seuilPval)), which(FC > 
        seuilFC))
    textseuil1 <- paste("Seuil FC = ", seuilFC, collapse = "", 
        sep = "")
    textseuil2 <- paste("p-value = ", seuilPval, collapse = "", 
        sep = "")
    UR <- paste("UR = ", length(droite), collapse = "", sep = "")
    DR <- paste("DR = ", length(gauche), collapse = "", sep = "")
    g <- rep(1, length(FC))
    g[gauche] <- 2
    g[droite] <- 3
    donnees <- data.frame(FC = FC, B = B, g = g)
    print(donnees[sample(1:length(FC), 10), ])
    graphes <- xyplot(B ~ FC, groups = g, panel = function(x, 
        y, subscripts, groups) {
        panel.abline(h = -log10(seuilPval), v = c(-seuilFC, seuilFC), 
            lty = 2)
        panel.superpose(x, y, groups = g, subscripts = subscripts, 
            col = c("black", "green", "red"))
    }, xlab = "log2(FC)", ylab = lab, main = titre, key = list(space = "top", 
        col = c("black", "black", "red", "green"), text = list(lab = c(textseuil1, 
            textseuil2, UR, DR)), columns = 2))
    plot(graphes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
