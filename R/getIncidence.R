##' Calculates the incidence of a gene set in an experiment 
##'
##' getIncidence returns an incidence vector of the location of the genes within a gene set in a list of genes in an experiment and vise-versa.
##' @title Determines the incidence of a gene set in a list of genes. 
##' @param universeIDs A vector containing the list of possible gene ids in the universe (or experiment).
##' @param set A GeneSet object containing a set of genes of interest
##' @return a list of inSet and inExp.  inSet is a vector with the same length as universeIDs.  Each value of inSet is 1 if the gene is in the set and 0 otherwise.  inExp is a vector with the same length as geneIds(set), the number of genes in the set.  Each value of inExp is 1 if the gene is in universeIDs and 0 otherwise.
##' @author Jessica L. Larson and Art Owen
##' @export
getIncidence <- function(universeIDs, set) {
    ##check if gene set
    if(is(set)[1] != "GeneSet" )
        {stop("set must be a GeneSet object") }     
    geneSetIds <- unlist(geneIds(set))
    inSet <- rep(0, length(universeIDs))
    inSet[match(geneSetIds, universeIDs)] <- 1
    inExp <- rep(0, length(geneSetIds))
    inExp[match(universeIDs, geneSetIds)] <- 1
    output <- list(inSet=inSet, inExp=inExp)
    output
}