##' Animation plot
##'
##' Animation plot
##' Generates plot in base R that gradually updates giving the impression of an animation
##' @title Animation plot
##' @param data Dimensionality reduction plot
##' @param ccm Dataframe of indices and momentums
##' @param delay Delay between frames in seconds
##' @param darken Passed to colGrad() function
##' @param lwd Line width
##' @param c.cex Size of poiints.
##' @param main Plot title
##' @param ... Passed to plot() function
##' @return Generates plot
##' @examples
##' \dontrun{
##' xx <- animPlot(x, ccm)
##' }
##' @author Wajid Jawaid
##' @importFrom graphics plot points lines
##' @export
animPlot <- function(data, ccm, delay = 1e-1, darken = 1, lwd = 1, c.cex = 1, main = "", ...) {
    plot(data, pch = 16, col = "gray", main = main, ...)
    ccm <- ccm[ccm[,1]!=0,]
    cc <- ccm[,1]
    ncols <- colGrad(length(cc), darken = darken)
    points(data[cc[1],,drop=FALSE], pch = 16, col = ncols[1], cex = c.cex)
    ## xr <- range(data[,1])
    ## yr <- range(data[,2])
    ## xp <- seq(xr[1], xr[2], length.out = 4)[1:2]
    ## yp <- seq(yr[1], yr[2], length.out = 4)[3:4]
    ## xin <- 0.02 * diff(xr)
    ## yin <- 0.02 * diff(yr)
    ## y2x <- xin / yin
    ## xp <- xp + xin
    ## yp <- yp + yin
    ## xm <- mean(xp)
    ## ym <- mean(yp)
    ## vec <- c(0,0)
    for (i in 2:length(cc)) {
        lines(data[c(cc[i-1],cc[i]),,drop = FALSE], col = ncols[i], lwd = lwd)
        points(data[cc[i],,drop=FALSE], pch = 16, col = ncols[i], cex = c.cex)
        ## if (any(vec != 0)) {
        ##     arrows(x0 = xm, y0 = ym, x1 = xm + vec[1], y1 = ym + vec[2], col = "white",
        ##            lwd = 2)
        ## }
        ## vec <- ccm[i, 2:3]
        ## vec[2] <- vec[2] / y2x
        ## vec <- unitLength(vec) / 10
        ## if (any(vec!=0)) {
        ##     arrows(x0 = xm, y0 = ym, x1 = xm + vec[1], y1 = ym + vec[2], col = "black",
        ##            lwd = 2)
        ##     }
        Sys.sleep(delay)
    }
}

##' Generates a GIF animation 
##'
##' Generates a GIF animation
##' @title Generates a GIF animating 
##' @param data Reduced dimensionality map to be used for visualisation
##' @param ccm Dataframe of indices and momentums
##' @param delay Delay between frames in seconds
##' @param darken Passed to colGrad() function
##' @param lwd Line width
##' @param c.cex Size of poiints.
##' @param main Title
##' @param gif Name of movie
##' @param img.name Name of temporary image files generated
##' @param plot.par Passed to R base par() function
##' @param point.col Colour of background points
##' @param arrowLength Modify length of arrow
##' @param ... Passed to plot() function
##' @return Produces an animated GIF with given file name
##' @examples
##' \dontrun{
##' xx <- animPlotGif(x, ccm, gif = "animation")
##' }
##' @author Wajid Jawaid
##' @importFrom animation saveGIF ani.options
##' @importFrom graphics segments arrows par
##' @export
animPlotGif <- function(data, ccm, delay = 1e-1, darken =1 , lwd =1, c.cex = 1, main = "", gif = "animation", img.name = "tempPlot", plot.par = NULL, point.col = "#333333", arrowLength = 0.1, ...) {
    ccm <- ccm[ccm[,1]!=0,]
    cc <- ccm[,1]
    ncols <- colGrad(length(cc), darken = darken)
    xr <- range(data[,1])
    yr <- range(data[,2])
    xp <- seq(xr[1], xr[2], length.out = 4)[1:2]
    yp <- seq(yr[1], yr[2], length.out = 4)[3:4]
    xin <- 0.02 * diff(xr)
    yin <- 0.02 * diff(yr)
    y2x <- xin / yin
    xp <- xp + xin
    yp <- yp + yin
    xm <- mean(xp)
    ym <- mean(yp)
    ani.options(autobrowse = FALSE, verbose = FALSE)
    ani.options(interval = delay)
    saveGIF(
        for (i in 1:length(cc)) {
            if (!is.null(plot.par)) par(plot.par)
            plot(data, pch = 16, col = point.col, main = main, ...)
            if (i > 1) {
                segments(data[cc[1:(i-1)],1], data[cc[1:(i-1)],2],
                         data[cc[2:i],1], data[cc[2:i],2], col = ncols[1:i], lwd = lwd)
            }
            points(data[cc[1:i],], pch = 16, col = ncols[1:i], cex = c.cex)
            ## theta <- atan2(ccm[i,3], ccm[i,2])
            vec <- ccm[i, 2:3]
            ## vec[2] <- vec[2] / y2x
            vec <- arrowLength * unitLength(vec)
            if (any(vec!=0)) {
                arrows(x0 = xm, y0 = ym, x1 = xm + vec[1], y1 = ym + vec[2], col = "white",
                       lwd = 2)
            }
        }, movie.name = paste0(gif, ".gif"), img.name = img.name, interval = delay)
}
