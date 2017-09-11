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
##' @param seed Set seed to aid reproducibility
##' @param drl_options Options passed to layout_with_drl()
##' @param filter.mean Default 1e-02. Threshold for mean filter
##' @param filter.cv Default 2. Threshold for CV/Fano filter
##' @param use_fano Default FALSE. If TRUE use Fano factor not CV.
##' @param verbose Default FALSE. If TRUE will return more details
##' @param oldOut Default NULL. Can return previous results for quicker analysis.
##' @param repeatDRL Default FALSE. If you wish to redo the layout using previous
##' iteration layout as starting point.
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
##' @importFrom stats dist prcomp runif
##' @importFrom rARPACK eigs
##' @export
goggles <- function(x, pcaDims = 90, nsig = 5, dmat = NULL, mkv = NULL, plotDims = 2,
                    kernSq = 2, seed = 0, drl_options = NULL, filter.mean = 0.01,
                    filter.cv = 2, use_fano = FALSE, verbose = FALSE, oldOut = NULL,
                    repeatDRL = FALSE, ...) {

    if (repeatDRL && is.null(oldOut)) {
        stop("Cannot repeat drL without previous iteration. Please supply in oldOut.")
    }
    
    if (!is.matrix(x)) {
        stop("x must be a matrix.")
    }

    sameSF <- FALSE
    pp <- oldOut$params
    pp.filter <- (pp$filter.mean != filter.mean ||
                  pp$filter.cv != filter.cv ||
                  pp$use_fano != use_fano)
    pp.filter <- !is.na(pp.filter) && pp.filter
    pp.pca <- pp$pcaDims != pcaDims
    pp.pca <- !is.na(pp.pca) && pp.pca
    if (is.na(oldOut$pca) && pp.pca) pp.pca <- FALSE
    
    if (is.null(dmat) && (!pp.filter && !pp.pca)) {
        dmat <- oldOut$dmat
        xpc <- oldOut$pca
        genesUsed <- oldOut$genesUsed
        sameSF <- TRUE                  # Same so Far
    } else if (is.null(dmat) && (pp.filter || pp.pca)) {
        cat("Preprocessing matrix ...")
        xx <- bgGeneNorm(x)
        genesUsed <- colnames(xx)
        xx <- filterGenes(xx, fano = use_fano, mu = filter.mean, cv = filter.cv, verbose = verbose)
        xx <- scale(xx)
        rm(xx); gc()
        cat("done.\nPerforming PCA ... ")
        xpc <- eigs(xx, which = "LR", k = pcaDims)$vectors
        cat("done.\nCalculating distance matrix ... ")
        dmat <- as.matrix(dist(xpc))
        cat("done.\n")
    } else {
        genesUsed <- NA
        xpc <- NA
    }

    pp.mkv <- !sameSF || (pp$kernSq != kernSq)
    
    if (is.null(mkv) && pp.mkv) {
        cat("Calculating variable sigmas ... ")
        sigmas <- calculateVariableSigmas(dmat, nsig) / kernSq
        cat("done.\nApplying gaussian kernel ... ")
        adj <- applyGaussianKernelwithVariableSigma(dmat, sigmas)
        cat("done.\n")
        diag(adj) <- 0
        adj <- adj / rowSums(adj)
    } else if (is.null(mkv) && !pp.mkv) {
        ajd <- oldOut$adj
        najd <- oldOut$nadj
        lmp <- oldOut$sparse[1]
        hmp <- oldOut$sparse[2]
        lvnClust <- oldOut$clust
        sameSF <- TRUE
    } else {
        adj <- mkv
    }
    

    ## summary(rowSums(adj>0))                 # Dense matrix?
    ## summary(colSums(adj>0))                 # Dense matrix?
    ## sum(colSums(adj)==0)                    # Any unreachable states
    ## summary(rowSums(adj))                   # Should be normalised

    if (!sameSF) {
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
        nadja <- nadj[allowableStates, allowableStates]
        cat("done.\nClustering ... ")
        lvnClust <- findLouvain(nadja)
    }

    sameSF <- sameSF && (pp$seed == seed)
    if (!sameSF && repeatDRL) {
        warning("Parameters have changed but previously calculated coordinates used as starting points. Was this intentional?")
    }
    
    if (sameSF && !repeatDRL) {
        l <- oldOut$l
    } else {
        cat("done.\nGenerating graph layout ... ")
        set.seed(seed)
        if (repeatDRL) {
            seedmat <- oldOut$l
        } else {
            seedmat <- matrix(runif(plotDims * nrow(nadja)), ncol = plotDims)
        }
        l <- layout_with_drl(lvnClust$gph, dim=plotDims, use.seed = 0, seed = seedmat,
                             options = drl_options)
        rownames(l) <- rownames(x)[allowableStates]
        colnames(l) <- paste0("drl_", 1:ncol(l))
    }
    cat("done.\n")

    return(list(l = l, clust = lvnClust, adj = adj, dmat = dmat,
                pca = xpc, sparse = c(lmp, hmp), nadj = nadj,
                nadja = nadja, seed = seed, genesUsed = genesUsed,
                params = list(pcaDims = pcaDims, nsig = nsig, kernSq = kernSq,
                              seed = seed, drl_options = drl_options,
                              filter.mean = filter.mean, filter.cv = filter.cv,
                              use_fano = use_fano)))
}
