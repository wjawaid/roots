##' Normalise by background gene set
##'
##' Normalise by background gene set.
##' Find background genes that are expressed at a lower percentage of
##' the total library size per cell than 'threshold' parameter. These
##' genes are used to calculate a normalisation factor.
##' @title Normalise by background gene set
##' @param x Matrix to be normalised with cells in rows and
##' genes in columns
##' @param threshold Default 0.05. The threshold below which a gene is
##' deemed background
##' @examples
##' \dontrun{
##' normGenes <- bgGeneNorm(x)
##' }
##' @return Returns a normalised matrix of same dimenions as 'x'
##' @author Wajid Jawaid
bgGeneNorm <- function(x, threshold = 0.05) {
    tc <- rowSums(x)
    tc <- tc %*% t(rep(1, ncol(x)))
    bgGenes <- apply(x < tc * threshold, 2, all)
    ftc <- rowSums(x[,bgGenes])
    x / ftc * mean(tc)    
}


##' Filter genes
##'
##' Filter genes
##' Filter genes by mean and either coefficient of variation, cv or
##' Fano factor.
##' @title Filter genes
##' @param mu Meam threshold
##' @param cv Coefficient of variation or Fano factor threshold. 
##' @param fano Default TRUE. Predicate treat CV as Fano factor or CV
##' @return Returns a filtered matrix with same number of cells but fewer
##' genes than 'x'
##' @examples
##' \dontrun{
##' expressionGenesFiltered <- filterGenes(x)
##' }
##' @author Wajid Jawaid
##' @inheritParams bgGeneNorm
##' @importFrom stats sd var
filterGenes <- function(x, mu = 0.01, cv = 2, fano = FALSE) {
    cmu <- colMeans(x)
    if (fano)
        ccv <- apply(x, 2, var) / cmu
    else
        ccv <- apply(x, 2, sd) / cmu
    plot(cmu, ccv)
    filteredGenes <- (cmu > mu) & (ccv > cv) 
    points(cmu[filteredGenes], ccv[filteredGenes], col=2)
    x <- x[,filteredGenes]
}
