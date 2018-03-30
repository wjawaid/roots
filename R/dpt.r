##' Calculate pseudotime
##'
##' Calculate diffusion pseudotime
##' @title Calculate pseudotime
##' @param dmap Output from diffuseMat function
##' @param startIndex Numerical index of starting cell
##' @return Vector of pseudotimes
##' @author Wajid Jawaid
##' @export
dpt <- function(dmap, startIndex) {
    vd <- (dmap$vectors - dmap$vectors[,startIndex])^2
    w <- dmap$values / (1 - dmap$values)
    pt <- t(t(vd) * w)
    rowSums(pt)
}
