##' View single cell dataset
##'
##' View single cell dataset
##' @title View single cell dataset
##' @param x Matrix with cells in rows and gene in columns
##' @param pcaDims Number of PCA dimensions to keep for distance measure
##' @param nsig Number of significant neighbours to keep for Gaussian kernel
##' @param dmat Optional. Give your own distance matrix
##' @param mkv Optional. Give your own markov matrix.
##' @param plotDims Default 2. Number of dimensions to plot
##' @param kernSq Factor to tighten kernel - operates on sigmas.
##' @param ... Additonal parameters not currently in use
##' @return A list of l, dimensionality reduced data.frame;
##' clust, returned from louvainClust();
##' adj, Sparse, pruned adjacency matrix;
##' dmat, distance matrix;
##' pca, PCA reduced matrix.
##' sparse, diagnostics on adj prior to applying sparseMarkov().
##' @examples
##' \dontrun{
##' xx <- goggles(x)
##' plot(xx$l)
##' }
##' @author Wajid Jawaid
##' @importFrom igraph layout_with_fr layout_with_drl
##' @importFrom stats dist prcomp
##' @export
goggles <- function(x, pcaDims = 90, nsig = 5, dmat = NULL, mkv = NULL, plotDims = 2,
                    kernSq = 2, ...) {
    if (!is.matrix(x)) {
        stop("x must be a matrix.")
    }
    
    if (is.null(dmat)) {
        cat("Preprocessind matrix ...")
        xx <- bgGeneNorm(x)
        xx <- filterGenes(xx, fano = FALSE)
        xx <- scale(xx)
        cat("done.\nPerforming PCA ... ")
        xpc <- prcomp(xx, center = FALSE, scale. = FALSE)$x[,1:pcaDims]
        cat("done.\nCalculating distance matrix ... ")
        dmat <- as.matrix(dist(xpc))
        cat("done.\n")
    }

    if (is.null(mkv)) {
        cat("Calculating variable sigmas ... ")
        sigmas <- calculateVariableSigmas(dmat, nsig) / kernSq
        cat("done.\nApplying gaussian kernel ... ")
        adj <- applyGaussianKernelwithVariableSigma(dmat, sigmas)
        cat("done.\n")
        diag(adj) <- 0
        adj <- adj / rowSums(adj)
    }
    ## summary(rowSums(adj>0))                 # Dense matrix?
    ## summary(colSums(adj>0))                 # Dense matrix?
    ## sum(colSums(adj)==0)                    # Any unreachable states
    ## summary(rowSums(adj))                   # Should be normalised

    lmp <- 1/min(apply(adj, 1, max))               # lowest max prob
    hmp <- 1/max(apply(adj, 1, max))               # highest max prob
    ## c(lmp, hmp)
    ## hist(1/apply(adj, 1, max), breaks= 100)
    cat("Making sparse ... ")
    adj <- sparseMarkov(adj, ceiling(lmp))
    ## summary(rowSums(adj>0))                 # More diagnostics
    ## summary(colSums(adj>0))
    ## sum(colSums(adj)==0)

    ## Take 2-step approach
    cat("done.\nPruning spurious edges ... ")
    nadj <- adj > 0
    storage.mode(nadj) <- "integer"
    nadj <- nadj %*% nadj
    nadj <- nadj >= 1
    nadj <- nadj * adj
    cat("done.\nExcluding non-allowable states ... ")
    
    ## Some states have no exit
    allowableStates <- rowSums(nadj) != 0
    sum(allowableStates)
    nadj <- nadj[allowableStates, allowableStates]
    cat("done.\nClustering ... ")
    
    lvnClust <- findLouvain(nadj)
    cat("done.\nGenerating graph layout ... ")
    l <- layout_with_drl(lvnClust$gph, dim=plotDims)
    rownames(l) <- rownames(x)[allowableStates]
    colnames(l) <- paste0("drl_", 1:ncol(l))
    cat("done.\n")
    return(list(l = l, clust = lvnClust, adj = adj, dmat = dmat,
                pca = xpc, sparse = c(lmp, hmp)))
}
