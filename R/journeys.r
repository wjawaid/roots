unitLength <- function(x) {
    vecLen <- sqrt(t(x) %*% x)
    if (vecLen != 0) return(x / vecLen)
    x
}

##' Find next cell function
##'
##' Find next cell function. Transitioin probabilities are modifed by calulating the cosine
##' of the angle between the current momentum vector and the vector on the rdmap required for
##' each transtion. The tranisiton probability is adjusted by multiplying by w1^(w2 * (cosine_angle)) and then normalising.
##' @title Find next cell function
##' @param rdmap reduced dimensionality matrix with cells in rows and dims in columns
##' @param tm Transition matrix
##' @param curInd Current state on tm
##' @param mom Current momentum vector
##' @param momAdj Weighting to adjust momentum. From 0-1. Lower numbers make smaller
##' adjustment to momentum vector.
##' @param w1 Parameter - Base used for modifying of tm probs.
##' @param w2 Parameter - Multiplifaction factor used for modifying tm probs.
##' @return Returns index of new cell and new momentum vector
##' @author Wajid Jawaid
fnc <- function(rdmap, tm, curInd, mom = NULL, momAdj = 0.5, w1 = exp(1), w2 = 1) {
    ## ca - cosine angle; rdmap - reduced dimensionality map, mom - momentum
    nCells <- nrow(rdmap)
    if (curInd == 0) return(sample(1:nCells, 1))
    if (is.null(mom) || all(mom == 0)) {
        cc <- sample(1:nCells, 1, prob = tm[curInd,])
        return(list(cc=cc, mom = rdmap[cc,] - rdmap[curInd,]))
    }
    ca <- t(t(rdmap) - rdmap[curInd,])
    ca <- apply(ca, 1, unitLength)
    ca <- as.numeric(mom %*% ca)
    ca <- w1^(w2*ca)
    ca <- tm[curInd,] * ca
    ca <- ca / sum(ca)
    cc <- sample(1:nCells, 1, prob = ca)
    newMom <- (rdmap[cc,] - rdmap[curInd,]) * momAdj
    newMom <- unitLength(newMom + (1 - momAdj) * mom)
    return(list(cc=cc, mom = newMom))
}

getTraj <- function(rdmap, tm, sourceCellInds, terminalCellsInd = NULL,
                    momAdj = 0.5, w1 = exp(1), w2 = 1, simLen = 50,
                    sim.seed = NULL) {
    if(!is.null(sim.seed)) set.seed(sim.seed)
    nCell <- vector("numeric", simLen)
    nCell[1] <- sample(sourceCellInds, 1)
    im <- matrix(0, ncol = ncol(rdmap), nrow = simLen)
    for (i in 2:(length(nCell))) {
        nc <- fnc(rdmap[,1:2], tm, nCell[i-1], im[i - 1,], momAdj = momAdj, w1 = w1, w2 = w2)
        nCell[i] <- nc[[1]]
        im[i,] <- nc[[2]]
        if (nCell[i] %in% terminalCellsInd) break
    }
    traj <- cbind(ind = nCell, mom.x = im[,1], mom.y = im[,2])
    attr(traj, "seed") <- sim.seed
    return(traj)
}
