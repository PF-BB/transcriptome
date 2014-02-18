afficheTable <-
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
        text(0.45/2, 0.9 - i * (0.9 - 0.01)/(n + 1), round((object)[i,1],2), 
            cex = cexval)
        polygon(c(0, 0.9), c(0.9 - (i + 0.5) * (0.9 - 0.01)/(n + 
            1), 0.9 - (i + 0.5) * (0.9 - 0.01)/(n + 1)))
        text(3 * 0.45/2, 0.9 - i * (0.9 - 0.01)/(n + 1), round((object)[i,2],2), 
            cex = cexval)
    }
    return(TRUE)
}
