##' Louvain clustering on transition matrix
##'
##' Louvain clustering on transition matrix
##' @title Louvain clustering on transition matrix
##' @param mkv Transition matrix
##' @return Returns a list with graph, dataframe and community object
##' @examples
##' \dontrun{
##' xx <- findLouvain(mkv)
##' xx$cll
##' }
##' @author Wajid Jawaid
##' @export
##' @importFrom igraph cluster_louvain communities graph.adjacency
findLouvain <- function(mkv) {
    gph <- graph.adjacency(mkv, weighted = TRUE, mode = "undirected")
    cll <- cluster_louvain(gph)
    lvnClust <- communities(cll)
    lvnClust <- lapply(1:length(lvnClust),
                       function(x) cbind.data.frame(cell=lvnClust[[x]],
                                                    class=as.character(x)))
    lvnClust <- do.call(rbind, lvnClust)
    rownames(lvnClust) <- lvnClust$cell
    lvnClust
    lvnClust[,1] <- NULL                # Remove duplicate name column
    lvnClust <- lvnClust[rownames(mkv),,drop=FALSE]
    return(list(gph=gph, clust=lvnClust, cll = cll))
    ## lvnClust <- lvnClust[rownames(pData(scd)),,drop=FALSE]    
}
