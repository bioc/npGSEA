 ##' Calculates an approximation of the permutation GSEA statistics and p-values
##'
##' This function calculates the permutation gene set enrichment analysis test statistic and p-value without actually running the permutation.  We account for the covariance among the genes within the set and approximate the corresponding permutation distribution.  For more details on the method see Larson and Owen (2014).
##' @title Runs the non-permutation GSEA 
##' @param x A matrix of expression data or an object of type ExpressionSet.  The columns of x represent samples in a given experiment.  The rows are genes.  The names of each row (or featureNames of the eSet) must be of the same type (e.g., entrez ids) as the ids of the gene set.
##' @param y A vector containing the treatment for each sample. The length of y must be more than 4 for the "chisq" approximation.
##' @param z A vector or matrix containing covariate(s) of interest, optional
##' @param set A GeneSet object containing a set of genes of interest or a GeneSetCollection object containing a collection of GeneSets
##' @param approx A string of either "norm" (default), "beta" or "chiSq".  If "norm", the normal approximation to the non-permutation GSEA is calculated and returned.  If "beta", the beta approximation is reported.  If "chiSq", the Chi-squared approximation to the permutation GSEA is calculated.
##' @param w A vector or list containing the weights of each gene in the set or sets, optional.  If w is a list, the number of elements in the list must correspond to the number of gene sets in the collection. 
##' @return an object with the corresponding GSEA results.  If approx="norm" an npGSEAResultNorm object is returned.  If approx="beta" a npGSEAResultBeta object is returned.  If approx="chiSq" a npGSEAResultChiSq object is returned.  If set is a GeneSetCollection (i.e., multiple sets of interest), then the corresponding npGSEAResultNormCollection, npGSEAResultBetaCollection, or npGSEAResultChiSqCollection is returned.
##' @author Jessica L. Larson and Art Owen
##' @export
npGSEA <- function(x, y, set, z = NULL, approx = c("norm", "beta", "chiSq"), w = NULL){
    approx <- match.arg(approx, c("norm", "beta", "chiSq"), several.ok = FALSE)    
    if(is(x)[1]== "ExpressionSet") { x <- exprs(x) }
    if(dim(x) [2] != length(y) )
        {stop("Must have the same amount of samples in x and y") }
        
    ##scale y to sum to zero
    y <- .adjustY(y)
      
    ######collection of sets:
    if (is(set)[1] == "GeneSetCollection"){
        if( (is.null(w) == FALSE) && ( is.list(w)== FALSE ) ){
            stop("You must provide weights for each gene set in 
            your collection as a list. 
            The length of w must be equal to the 
            length of your GeneSetCollection")
        }
        if( (is.null(w) == FALSE) && ( length(w)!=length(set) ) ){
            stop("You must provide weights for each gene 
            set in your collection.  
            The length of w is not equal to the 
            length of your GeneSetCollection")
        }
        ##make an empty list of weights if w is null
        if (is.null(w) ==TRUE )  { w <- vector("list", length(set)) }
        ##norm approx
        if (approx=="norm"){
            res <- mapply(function(singleSet, singleW){
            ##prep data
            xyz <- .prepXYZ(x, y,  z, singleSet)
            xg <- xyz$xg
            y <- xyz$y
            wg <- .prepW(singleW, singleSet,  xyz$inset)
            ##run analysis
            runNormApprox(xg, y, wg, singleSet)            
            }, set, w )
            output <- new("npGSEAResultNormCollection", res)
        }
        ##beta approx
        if (approx=="beta"){              
            res <- mapply(function(singleSet, singleW){
            xyz <- .prepXYZ(x, y,  z, singleSet)
            xg <- xyz$xg
            y <- xyz$y
            wg <- .prepW(singleW, singleSet,  xyz$inset)
            runBetaApprox(xg, y, wg, singleSet)            
            }, set, w )
            output <- new("npGSEAResultBetaCollection", res)
        }
        ##chisq approx
        if (approx=="chiSq"){
            res <- mapply(function(singleSet, singleW){
            xyz <- .prepXYZ(x, y,  z, singleSet)
            xg <- xyz$xg
            y <- xyz$y
            wg <- .prepW(singleW, singleSet, xyz$inset)
            runChisqApprox(xg, y, wg, singleSet)            
            }, set, w )
            output <- new("npGSEAResultChiSqCollection", res)
        }
    }  	

    #########single gene sets:
    else if (is(set)[1] == "GeneSet"){
        ##prep data
        xyz <- .prepXYZ(x, y,  z, set)
        xg <- xyz$xg
        y <- xyz$y
        wg <- .prepW(w, set, xyz$inset)
        ###run the analysis:
        if (approx=="norm"){
            output <- runNormApprox(xg, y, wg, set) 
        }	
        if (approx=="beta"){
            output <- runBetaApprox(xg, y, wg, set) 
        }	
        if (approx=="chiSq"){
            output <- runChisqApprox(xg, y, wg, set) 
        }
    }
    output
}

