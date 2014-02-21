plotQC2 <-
function (x, fc.line.col = "black", sf.ok.region = "light blue", 
    chip.label.col = "black", sf.thresh = 3, gdh.thresh = 1.58, 
    ba.thresh = 2.32, present.thresh = 10, bg.thresh = 20, label = NULL, 
    main = "QC Stats", usemid = F, spread = c(-8, 8), type = "l", 
    cex = 1, ...) 
{
    if (type == "c") {
        .plot.qc.stats2(x, fc.line.col, sf.ok.region, chip.label.col, 
            sf.thresh, gdh.thresh, ba.thresh, present.thresh, 
            bg.thresh, label, main, usemid, cex, ...)
        return()
    }
    old.par <- par()
    par(mai = c(0, 0, 0, 0))
    sfs <- log2(sfs(x))
    n <- length(sfs)
    meansf <- mean(sfs)
    dpv <- percent.present(x)
    dpv <- (round(100 * dpv))/100
    abg <- avbg(x)
    abg <- (round(100 * abg))/100
    if (is.null(label)) {
        label <- names(maxbg(x))
    }
    d1 <- 0
    d2 <- 0
    d3 <- 0
    for (i in 1:n) {
        for (j in 1:n) {
            d1 <- max(abs(sfs[i] - sfs[j]), d1)
            d2 <- max(abs(dpv[i] - dpv[j]), d2)
            d3 <- max(abs(abg[i] - abg[j]), d3)
        }
    }
    m <- matrix(c(4, 2, 1, 3), nrow = 2, ncol = 2)
    layout(m, c(1, 2), c(0.1, 1))
    if (is.null(main)) {
        main = ""
    }
    plot(0, 0, xlim = range(0, 1), ylim = range(0, 1), type = "n", 
        yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
    text(0.5, 0.5, labels = main, adj = 0, cex = cex * 2)
    plot(0, 0, xlim = range(0, 1), ylim = range(-1, n), type = "n", 
        yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
    text(1, (1:n) - 0.5, labels = label, adj = 1, cex = cex)
    plot(0, 0, xlim = spread, ylim = c(-1, n), type = "n", xaxs = "i", 
        yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
    x1 <- (sf.thresh/2 + meansf)
    y1 <- 0
    x2 <- (-sf.thresh/2 + meansf)
    y2 <- n
    polygon(c(x1, x2, x2, x1), c(y1, y1, y2, y2), col = sf.ok.region, 
        border = sf.ok.region)#
    lines(c(0, 0), c(0, n), lty = 1, col = fc.line.col)
    lines(c(-1, -1), c(0, n), lty = 2, col = "grey")
    lines(c(-2, -2), c(0, n), lty = 2, col = "grey")
    lines(c(-3, -3), c(0, n), lty = 2, col = fc.line.col)
    lines(c(1, 1), c(0, n), lty = 2, col = "grey")
    lines(c(2, 2), c(0, n), lty = 2, col = "grey")
    lines(c(3, 3), c(0, n), lty = 2, col = fc.line.col)
    lines(c(1.58, 1.58), c(0, n), lty = 2, col = 2)
    lines(c(2.32, 2.32), c(0, n), lty = 1, col = 2)
    text(3, -1, "3", pos = 3, col = fc.line.col, cex = cex)
    text(2, -1, "2", pos = 3, col = fc.line.col, cex = cex)
    text(1, -1, "1", pos = 3, col = fc.line.col, cex = cex)
    text(-3, -1, "-3", pos = 3, col = fc.line.col, cex = cex)
    text(-2, -1, "-2", pos = 3, col = fc.line.col, cex = cex)
    text(-1, -1, "-1", pos = 3, col = fc.line.col, cex = cex)
    text(0, -1, "0", pos = 3, col = fc.line.col, cex = cex)
    text(-8, -0.6,"PP : % presents", pos = 4, cex = cex)
    text(-8, -0.8,"AvgBG : Moyenne Background", pos = 4, cex = cex)
    rats <- ratios(x)
    if (!usemid) {
        gdh <- rats[, 2]
        ba <- rats[, 1]
    }
    else {
        gdh <- rats[, 4]
        ba <- rats[, 3]
    }
    bb <- x@bioBCalls
    for (i in 1:n) {
        x1 <- spread[1]
        x2 <- spread[2]
        y1 <- i - 1
        y2 <- i - 1
        lines(c(x1, x2), c(y1, y2), lty = 2, col = "light grey")
        if (d1 > sf.thresh) {
            col = "red"
        }
        else {
            col = "blue"
        }
        x1 <- sfs[i]
        y1 <- i - 0.25
        lines(c(0, x1), c(y1, y1), col = col)
        points(x1, y1, col = col, pch = 20)
        x2 <- gdh[i]
        y2 <- i - 0.5
        if (gdh[i] > gdh.thresh) {
            col = "red"
        }
        else {
            col = "blue"
        }
        points(x2, y2, pch = 1, col = col)
        x2 <- ba[i]
        y2 <- i - 0.5
        if (ba[i] > ba.thresh) {
            col = "red"
        }
        else {
            col = "blue"
        }
        points(x2, y2, pch = 2, col = col)
        if (d2 > present.thresh) {
            col = "red"
        }
        else {
            col = "blue"
        }
        x2 <- spread[1]
        y2 <- i - 0.25
        dpvs <- paste(dpv[i], "% (PP)", sep = "")
        text(x2, y2, label = dpvs, col = col, pos = 4, cex = cex)
        if (d3 > bg.thresh) {
            col = "red"
        }
        else {
            col = "blue"
        }
        x2 <- spread[1]
        y2 <- i - 0.75
        abgs <- paste(abg[i], " (AvgBG)", sep = "")
        text(x2, y2, label = abgs, col = col, pos = 4, cex = cex)
        if (bb[i] != "P") {
            text(0, i - 1, label = "bioB", col = "red", cex = cex)
        }
    }
    plot(0, 0, xlim = range(0, 1), ylim = range(0, 1), type = "n", 
        yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
    if (!usemid) {
        points(0.25, 0.2, pch = 1)
        text(0.3, 0.2, colnames(rats)[2], pos = 4, cex = cex)
        points(0.25, 0.4, pch = 2)
        text(0.3, 0.4, colnames(rats)[1], pos = 4, cex = cex)
        points(0.25, 0.6, pch = 20)
        text(0.3, 0.6,"Scale Factor", pos = 4, cex = cex)
    }
    else {
        points(0.25, 0.2, pch = 1)
        text(0.3, 0.2, colnames(rats)[4], pos = 4, cex = cex)
        points(0.25, 0.4, pch = 2)
        text(0.3, 0.4, colnames(rats)[3], pos = 4, cex = cex)
        points(0.25, 0.6, pch = 20)
        text(0.3, 0.6,"Scale Factor", pos = 4, cex = cex)
    }
    ow <- options("warn")$warn
    options(warn = -1)
    par(old.par)
    options(warn = ow)
}
