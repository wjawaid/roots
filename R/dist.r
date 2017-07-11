##' Fast vectorised Euclidean distance calculator
##'
##' Calculates Euclidean distances between vectors arranged as columns in a matrix.
##' @title Fast vectorised Euclidean distance calculator
##' @param x Matrix with vectors in columns.
##' @param squared Will not perform the square root, i.e. will return the squared `L2-norm'.
##' @examples
##' \dontrun{
##' dist <- fastDist(x)
##' }
##' @return Returns a matrix of pairwise distances
##' @author Wajid Jawaid
##' @export
fastDist <- function(x, squared = FALSE) {
    a <- colSums(x^2)
    a <- a * matrix(1, ncol(x), ncol(x))
    a <- a + t(a)
    ab <- t(x) %*% x
    d <- a - 2 * ab
    diag(d) <- 0
    if (!squared) d <- sqrt(d)
    return(d)
}

##' Make markov matrix sparse
##'
##' Make markov matrix sparse
##' Choose knn as the maximum number of similar cells are likely to exist in your dataset.
##' @title Make markov matrix sparse
##' @param mkv Markov matric
##' @param knn Number of nearest neighbours. See above.
##' @return Markovian sparse matrix.
##' @author Wajid Jawaid
sparseMarkov <- function(mkv, knn) {
    predicateFilter <- mkv <= (1/knn)
    mkv[predicateFilter] <- 0
    rSums <- rowSums(mkv) == 0
    if (any(rSums)) {
        cat(sum(rSums), " cells have no transitions. Try increasing knn.\n")
    }
    mkv <- mkv / rowSums(mkv)
}
