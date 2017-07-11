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
##' @param varEst Number of alternatives to sample for estimating variance.
##' @examples
##' \dontrun{
##' nextCell <- fnc(rdmap, tm, curInd)
##' }
##' @return Returns index of new cell and new momentum vector
##' @author Wajid Jawaid
##' @export
fnc <- function(rdmap, tm, curInd, mom = NULL, momAdj = 0.5, w1 = exp(1), w2 = 1, varEst = 10) {
    ## ca - cosine angle; rdmap - reduced dimensionality map, mom - momentum
    nCells <- nrow(rdmap)
    if (curInd == 0) return(sample(1:nCells, 1))
    if (is.null(mom) || all(mom == 0)) {
        cc <- sample(1:nCells, varEst + 1, prob = tm[curInd,], replace = TRUE)
        return(list(cc=cc[1], mom = rdmap[cc[1],] - rdmap[curInd,], altCells = cc[-1]))
    }
    ca <- t(t(rdmap) - rdmap[curInd,])
    if (any(is.na(ca))) {
        errFile <- file("Error.log")
        writeLines(paste0("1. NA introduced on calculating transtion vectors.\n",
                          "Current Index: ", curInd,
                          "\nMomentum", mom, "\n"),
                   con = errFile)
        flush(errFile)
        close(errFile)
    }
    ## if (any(is.na(ca))) writeLines(paste0("1. NA introduced on calculating transtion vectors.\n", "Current Index: ", curInd, "\nMomentum", mom, "\n"), file = "Error.log")
    ca <- apply(ca, 1, unitLength)
    if (any(is.na(ca))) {
        errFile <- file("Error.log")
        writeLines(paste0("2. NA introduced on making unit vectors.\n",
                          "Current Index: ", curInd,
                          "\nMomentum", mom, "\n"),
                   con = errFile)
        flush(errFile)
        close(errFile)
    }
    ## if (any(is.na(ca))) writeLines(paste0("2. NA introduced on making unit vectors.\n", "Current Index: ", curInd, "\nMomentum", mom, "\n"), file = "Error.log")
    ca <- as.numeric(mom %*% ca)
    if (any(is.na(ca))) {
        errFile <- file("Error.log")
        writeLines(paste0("3. NA introduced on calculating dot product.\n", 
                          "Current Index: ", curInd,
                          "\nMomentum", mom, "\n"),
                   con = errFile)
        flush(errFile)
        close(errFile)
    }
    ## if (any(is.na(ca))) writeLines(paste0("3. NA introduced on calculating dot product.\n", "Current Index: ", curInd, "\nMomentum", mom, "\n"), file = "Error.log")
    ca <- w1^(w2*ca)
    if (any(is.na(ca))) {
        errFile <- file("Error.log")
        writeLines(paste0("4. NA introduced on converting angle to weighting factor\n",
                          "Current Index: ", curInd,
                          "\nMomentum", mom, "\n"),
                   con = errFile)
        flush(errFile)
        close(errFile)
    }
    
    ca2 <- tm[curInd,] * ca
    ## if (any(is.na(ca2))) {
    ##     cat("5\n")
    ##     cat(curInd, "\n")
    ##     cat(any(is.na(tm[curInd,])))
    ##     cat(tm[curInd,], "\n")
    ##     cat(ca, "\n")
    ## }
    ca <- ca2 / sum(ca2)
    ## if (any(is.na(ca))) Error(paste0("5. NA introduced on dividing\n", "Current Index: ", curInd, "\nMomentum", mom, "\nSumCA2", sum(ca2), "\n"))
    cc <- sample(1:nCells, varEst + 1, prob = ca, replace = TRUE)
    ## tryCatch(
    ##     cc <- sample(1:nCells, varEst + 1, prob = ca, replace = TRUE)
    ## , error = function(e) {
    ##     errFile <- file("Error.log")
    ##     writeLines(paste0("Current Index: ", curInd, "\nMomentum: ", paste(mom, collapse = ","), "\nSumCA2: ", sum(ca2)), con = errFile)
    ##     flush(errFile)
    ##     close(errFile)
    ## })
    
    newMom <- (rdmap[cc[1],] - rdmap[curInd,]) * momAdj
    newMom <- unitLength(newMom + (1 - momAdj) * mom)
    return(list(cc=cc[1], mom = newMom, altCells = cc[-1]))
}





##' Return a plausible developmental journey
##'
##' Return a plausible developmental journey
##' @title Find a plausible developmental journey
##' @param sourceCellInds Starting sell indices
##' @param terminalCellsInd Terminal cell indices
##' @param simLen Maximum number of allowable tranisitons
##' @param sim.seed Random seed for reproducibility
##' @return Returns a data.frame of ordered indices and momentums
##' @examples
##' \dontrun{
##' traj <- getTraj(rdmap, tm, startCells, terminalCells)
##' }
##' @inheritParams fnc
##' @author Wajid Jawaid
##' @export
getTraj <- function(rdmap, tm, sourceCellInds, terminalCellsInd = NULL,
                    momAdj = 0.5, w1 = exp(1), w2 = 1, simLen = 50,
                    sim.seed = NULL, varEst = 10) {
    if(!is.null(sim.seed)) set.seed(sim.seed)
    nCell <- vector("numeric", simLen)
    nCell[1] <- sample(sourceCellInds, 1)
    im <- matrix(0, ncol = ncol(rdmap), nrow = simLen)
    altCells <- matrix(0, ncol = varEst, nrow = simLen)
    colnames(altCells) <- paste0("alt", 1:varEst)
    for (i in 2:(length(nCell))) {
        nc <- fnc(rdmap[,1:2], tm, nCell[i-1], im[i - 1,], momAdj = momAdj, w1 = w1, w2 = w2, varEst = varEst)
        nCell[i] <- nc[[1]]
        im[i,] <- nc[[2]]
        altCells[i,] <- nc[[3]]
        if (nCell[i] %in% terminalCellsInd) break
    }
    traj <- cbind(ind = nCell, mom.x = im[,1], mom.y = im[,2], altCells)
    attr(traj, "seed") <- sim.seed
    return(traj)
}
