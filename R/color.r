##' Generates a smooth colour gradient
##'
##' Generates a smooth colour gradient
##' Goes from red to red/green to green to green/blue to blue to blu/red 
##' @title Generates a smooth colour gradient
##' @param x Number of colours required
##' @param darken Multiplication factor. Must be less than 1. Smaller the darker.
##' @return Returns vector of RGB colours
##' @examples
##' gradientColors <- colGrad(10)
##' @author Wajid Jawaid
##' @importFrom grDevices rgb
##' @export
colGrad <- function(x, darken = 1) {
    clts <- data.frame(R = c(255, 255, 0,   0,   0,   255),
                       G = c(0,   255, 255, 255, 0,   0),
                       B = c(0,   0,   0,   255, 255, 255))
    if (x <= 6) {
        rcols <- clts[1:x,]
        rcols <- rcols / 255
        return(apply(rcols, 1, function(x) rgb(x[1], x[2], x[3])))
    }
    cls <- data.frame(R = rep(0, x), G = 0, B = 0)
    cls[1,1] <- 255
    cj <- (x - 6) %/% 5
    cj <- rep(cj, 5)
    cr <- (x - 6) %% 5
    if (cr != 0) cj[1:cr] <- cj[1:cr] + 1
    if (sum(cj) + 6 != x) stop("Error 1.")
    bks <- c(1, 1+cumsum(cj + 1))
    cls[bks,] <- clts
    for (i in 1:5) {
        if (cj[i] !=0) {
            nc <- bks[i + 1] - bks[i]
            dt <- cls[c(bks[i],bks[i+1]), ]
            cc <- 0:(nc-1) %*%
                matrix(apply(dt, 2, function(x) (x[2] - x[1]) / nc), 1)
            cc <- t(t(cc) + as.numeric(dt[1,]))
            cls[(bks[i] + 1):(bks[i+1] - 1), ] <- cc[-1,]
        }
    }
    rcols <- cls / 255
    rcols[rcols>1] <- 1
    rcols <- rcols * darken
    return(apply(rcols, 1, function(x) rgb(x[1], x[2], x[3])))
}
