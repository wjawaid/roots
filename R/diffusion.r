##' Generic diffusion function using automated
##' individualised sigma calculation
##'
##' Generic diffusion function using automated individualised sigma calculation.
##'
##' A Gaussian kernel is applied to the chosen distance metric producing
##' an \eqn{n \times n} square unnormalised symmetric transition matrix, \eqn{A}.
##' Let \eqn{D} be an \eqn{n \times n} diagonal matrix with row(column) sums of
##' \eqn{A} as entries. The density corrected transition matrix will now
##' be:
##' 
##' \deqn{D^{-1} A D^{-1}}{D^{-1} * A * D^{-1}}
##'
##' and can be normalised:
##' 
##' \deqn{B^{-1} D^{-1} A D^{-1}}{B^{-1} * D^{-1} * A * D{^-1}}
##'
##' where \eqn{B} is an \eqn{n \times n} diagonal matrix with row sums of
##' the density corrected transition matrix as entries. The eigen decomposition of
##' this matrix can be simplified by solving the symmetric system:
##'
##' \deqn{B^{-\frac{1}{2}} D^{-1} A D^{-1} B^{-\frac{1}{2}} R^\prime = %
##'       R^\prime \lambda^\prime}{%
##'       B^{-0.5} * D^{-1} * A * D^{-1} * B^{-0.5} = R^\prime * L^\prime}
##'
##' where \eqn{R^\prime}{R^\prime} is a matrix of the right eigenvectors that solve
##' the system and \eqn{\lambda^\prime}{L'} is the corresponding eigenvalue
##' diagonal matrix. Now the solution of:
##'
##' \deqn{B^{-1} D^{-1} A D^{-1} R = R \lambda}{%
##'       B^{-1} * D^{-1} * A * D^{-1} * R = R * L}
##'
##' in terms of \eqn{R^\prime} and \eqn{B^{-\frac{1}{2}}} is:
##'
##' \deqn{B^{-1} D^{-1} A D^{-1} B^{-\frac{1}{2}} R^\prime = %
##'       B^{-\frac{1}{2}} R^\prime \lambda^\prime}{%
##'       B^{-1} * D^{-1} * A * D^{-1} * B^{-0.5} = B^{-0.5} * R' * L'}
##'
##' and
##'
##' \deqn{R = B^{-\frac{1}{2}} R^\prime}{R = B^{-0.5} * R'}
##'
##' This \eqn{R} without the first eigen vector is returned as the diffusion map.
##' @title Generic diffusion function
##' @param data Matrix of data with genes in rows and cells in columns.
##' @param ndims Number of dimensions to return
##' @param nsig For automatic sigma calculation
##' @param removeFirst Default TRUE. Removes the first eigenvector
##' @param useARPACK Default TRUE. Uses Arnoldi algorithm for eignvector calculations
##' @param distfun A different distance function that returns the \strong{squared}
##' distance
##' @param sigmas Manually provide sigma
##' @param sqdistmat \strong{Squared} distance matrix.
##' Give your own squared distance matrix.
##' @return List output containing:
##' \tabular{rl}{%
##'   \emph{values} \tab Eigenvalues, excluding the first eigenvalue, which should%
##'                       always be 1.\cr
##'    \emph{vectors} \tab Matrix of eigen vectors in columns, first eigen vector%
##'                        removed.\cr
##'    \emph{nconv} \tab Number of eigen vectors/values that converged.\cr
##'    \emph{niter} \tab Iterations taken for Arnoldi algorithm to converge.\cr
##'    \emph{nops} \tab  Number of operations. \cr
##'    \emph{val0} \tab 1st eigen value - should be 1. If not be suspicious!\cr
##'    \emph{vec0} \tab 1st eigen vector - should be \eqn{n^{-\frac{1}{2}}}{1/sqrt(n)},%
##'                      where n is the number of cells/samples.\cr
##'    \emph{usedARPACK} \tab Predicates use of ARPACK for spectral decomposition.\cr
##'    \emph{distfun} \tab Function used to calculate the squared distance.\cr
##'    \emph{nn} \tab Number of nearest neighbours used for calculating \code{sigmas}.\cr
##'    \emph{d2} \tab Matrix of squared distances, returned from \code{distfun}.\cr
##'    \emph{sigmas} \tab Vector of sigmas. Same length as number of cells if individual\cr
##'                  \tab sigmas were calculated, otherwise a scalar if was supplied.\cr
##'    \emph{gaussian} \tab Unnormalised transition matrix after applying Gaussian.\cr
##'    \emph{markov}  \tab Normalised \code{gaussian} matrix.\cr
##'    \emph{densityCorrected} \tab Matrix after applying density correction to %
##'                                 \code{markov}.\cr
##' }
##' @examples
##' \dontrun{
##' xx <- diffuseMat(x)
##' }
##' @author Wajid Jawaid
##' @importFrom rARPACK eigs
##' @references
##' Haghverdi, L., Buettner, F., Theis, F.J., 2015. Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics 31, 2989–2998.
##'
##' Haghverdi, L., Büttner, M., Wolf, F.A., Buettner, F., Theis, F.J., 2016. Diffusion pseudotime robustly reconstructs lineage branching. Nat Meth 13, 845–848.
##'
##' Angerer, P., Haghverdi, L., Büttner, M., Theis, F.J., Marr, C., Buettner, F., 2016. destiny: diffusion maps for large-scale single-cell data in R. Bioinformatics 32, 1241–1243.
##' @export
diffuseMat <- function(data, ndims = 20, nsig = 5,
                        removeFirst = TRUE, useARPACK = TRUE,
                        distfun = NULL, sigmas = NULL, sqdistmat = NULL) {
    if (!is.null(sqdistmat)) {
        d2 <- sqdistmat
    } else {
        cat("Calculating distance matrix. Please wait... ")
        if (is.null(distfun)) {
            d2 <- as.matrix(fastDist(data, squared = TRUE))       # distance matrix calculation
        } else d2 <- as.matrix(distfun(data))
        cat("Done.\n")
    }
    if (is.null(sigmas)) {
        cat("Calculating sigmas. Please wait... ")
        sigmas <- apply(d2, 1, function(x) sqrt(sort(x)[nsig])/2)
        cat("Done.\n")
        cat("Applying Gaussian kernel.\n")
        W <- applyGaussianKernelwithVariableSigma(d2, sigmas)
    } else {
        W <- exp(-d2 / (2*sigmas^2))
    }
    diag(W) <- 0
    cat("Calculating normalised transition matrix.\n")
    markov <- W / rowSums(W)
    cat("Performing density correction.\n")
    D <- rowSums(W)
    q <- D %*% t(D)                     # Calculate local density
    H <- W / q                          # Correct for local density
    dH <- rowSums(H)
    ## cat("Calculating density corrected normalised transition matrix.\n")
    ## markovDensityCorrected <- H / dH
    cat("Calculating related symmetric system.\n")
    rdH <- 1/sqrt(dH)
    Hp <- H * (rdH %*% t(rdH))
    cat("Calculating eigen vectors for symmetric system. Please wait... ")
    n <- nrow(d2)
    # Sub-function to sort and select largest eigenvecs/vals
    if (useARPACK) {                    # Eigen decomposition
        decomp <- rARPACK::eigs(Hp, which = "LR", k = ndims + 1)
    } else {   
        decomp <- eigen(Hp)
    }
    colnames(decomp$vectors) <- paste0("DC", 0:(ncol(decomp$vectors)-1))
    cat(" Done.\nTransforming to eigen vectors of normalised transition matrix ... ")
    decomp$vectors <- decomp$vectors * rdH
    rownames(decomp$vectors) <- colnames(data)
    ## decomp$vectors <- decomp$vectors / sqrt(colSums(H))
    cat(" Done.\nChoosing dims... ")

    if (removeFirst) {
        startDim <- 2
        decomp$val0 <- Re(decomp$values[1])
        decomp$vec0 <- Re(decomp$vectors[,1])
        decomp$values <- Re(decomp$values[-1])
        decomp$vectors <- Re(decomp$vectors[,-1])
    } else startDim <- 1
    cat(" Done.\n")
    if (!useARPACK) decomp$nconv <- decomp$niter <- integer(0)
    decomp$usedARPACK <- useARPACK
    decomp$distfun <- distfun
    decomp$nn <- nsig
    decomp$d2 <- d2
    decomp$sigmas <- sigmas
    rownames(W) <- colnames(W)
    decomp$gaussian <- W
    rownames(markov) <- colnames(markov)
    decomp$markov <- markov
    rownames(H) <- colnames(H)
    decomp$densityCorrected <- H
    return(decomp)
}

##' Predicts diffusion map projection from new data points
##'
##' Predicts diffusion map projection from new data points
##' @title Predicts diffusion map projection from new data points
##' @param dm Output from diffuseMat2 function
##' @param x Matrix of new data points. Features in rows and cells in
##' columns.
##' @param data Original data used to generate diffusion map
##' @param distfun A distance function that takes new data as first
##' paramter and previous data as second variable returning a squared
##' distance measure, with each sample in the rows and distance to
##' previous data points in columns, e.g. function(x, y) (1 - cor(x, y))^2.
##' @return Returns a matrix with projected diffusion components.
##' @examples
##' \dontrun{
##' y <- diffuseProj(xx, newData, oldData, function(z) (1-cor(z))^2)
##' }
##' @author Wajid Jawaid
##' @export
diffuseProj <- function(dm, x, data, distfun) {
    ## Calculate new distances
    d2 <- distfun(x, data)             # dim(d2) is ncols(x) by ncols(data)
    if (length(dm$sigmas) == 1) {
        W <- exp(-d2 / (2*dm$sigmas^2))    # Gaussian Kernel
    } else {
        #Find distance to nsig cell then / 2.
        sigmas <- apply(d2, 1, function(x) sqrt(sort(x)[dm$nn])/2)
        W <- applyGaussianKernelwithVariableSigma(d2, sigmas, dm$sigmas)
    }
    ## diag(W) <- 0                       # No self-transitions
    W[which(W == 1)] <- 0
    ## Perform density normalisation by dividing by the product of the sums of
    ## incoming and outgoing weights
    H <- W / rowSums(W)
    H <- t(t(H) / rowSums(dm$gaussian))
    ## Make markovian
    Hp <-  H / rowSums(H)
    predMat <- Hp %*% cbind(dm$vec0, dm$vectors) %*% diag(c(1, 1/dm$values))
    colnames(predMat) <- paste0("DC", 0:(ncol(predMat)-1))
    rownames(predMat) <- colnames(x)
    predMat[,-1, drop=FALSE]
}

##' Calculates sigmas for a distance matrix
##'
##' Calculates sigmas for a distance matrix
##' Using Laleh Hagherverdi's method
##' @title Calculates sigmas for a distance matrix
##' @param d Square distance matrix with 0 diagonal
##' @param knn Number of nearest neighbours to use for calculation
##' @examples
##' \dontrun{
##' sigmas <- calculateVariableSigmas(dist, 5)
##' }
##' @return Returns a vector of sigmas
##' @author wj241
calculateVariableSigmas <- function(d, knn) {
    if (nrow(d) != ncol(d)) stop("Distance matrix not square.")
    if (all(diag(d) != 0)) stop("Self loops present")
    return(apply(d, 1, function(x) sqrt(sort(x)[knn])/2))
}

##' Apply Gaussian Kernel using Laleh Haghverdi's variable sigma
##'
##' Apply Gaussian Kernel using Laleh Haghverdi's variable sigma
##' @title Apply Gaussian Kernel using Laleh Haghverdi's variable sigma
##' @param d2 Squared distance metric
##' @param rsigmas Sigmas for cells in the rows
##' @param csigmas Sigmas for cells in the columns
##' @return Returns matrix of same size as d2.
##' @examples
##' \dontrun{
##' d <- applyGaussianKernelwithVariableSigma(dist, sigmas)
##' }
##' @author Wajid Jawaid
applyGaussianKernelwithVariableSigma <- function(d2, rsigmas, csigmas = NULL) {
    if (is.null(csigmas)) {
        sigmaMat <- rsigmas %*% t(rsigmas)
        sigmaSqd <- rsigmas^2
        sigmaSum <- expand.grid(1:length(sigmaSqd), 1:length(sigmaSqd))
        sigmaSum <- sigmaSqd[sigmaSum[,1]] + sigmaSqd[sigmaSum[,2]]
        sigmaSum <- matrix(sigmaSum, length(sigmaSqd))
    } else {
        sigmaMat <- rsigmas %*% t(csigmas)
        rsigmaSqd <- rsigmas^2
        csigmaSqd <- csigmas^2
        sigmaSum <- expand.grid(1:length(rsigmaSqd), 1:length(csigmaSqd))
        sigmaSum <- rsigmaSqd[sigmaSum[,1]] + csigmaSqd[sigmaSum[,2]]
        sigmaSum <- matrix(sigmaSum, length(rsigmaSqd))
    }
    W <- sqrt(2*sigmaMat / sigmaSum) * exp(-d2 / sigmaSum)
}
