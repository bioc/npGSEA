\name{npGSEA}
\alias{npGSEA}
\title{Calculates an approximation of the permutation GSEA statistics and p-values}
\usage{
  npGSEA(x, y, set, covars = NULL, approx = c("norm", "beta", "chiSq"), w = NULL, epsilonBetaAdj=TRUE, scaleXY=TRUE )
}
\arguments{
  \item{x}{A matrix of expression data or an object of type ExpressionSet.  
      The columns of x represent samples in a given experiment.  The rows are genes.  
      The names of each row (or \code{featureNames} of the eSet) 
      must be of the same type (e.g., entrez ids) as the ids of the gene set.}

  \item{y}{A vector containing the treatment for each sample. 
      The length of y must be more than 4 for the "chisq" approximation.  
      Each treatment group must have at least two observations for all approximation methods.  
      There can only be two treatment groups.}
  
  \item{covars}{A vector or matrix containing covariate(s) of interest, optional}
  
  \item{set}{A \code{GeneSet} object containing a set of genes of interest or a 
  	\code{GeneSetCollection} object containing a collection of GeneSets}

  \item{approx}{A string of either "norm" (default), "beta" or "chiSq".  
      If "norm", the normal approximation to the non-permutation GSEA is calculated and returned.  
      If "beta", the beta approximation is reported.  If "chiSq", the Chi-squared approximation to the 
      permutation GSEA is calculated.}

  \item{w}{A vector or list containing the weights of each gene in the set or sets, optional.  If w is a list, 
  	the number of elements in the list must correspond to the number of gene sets in the collection.}
  	
  	  \item{epsilonBetaAdj}{A boolean indicating whether or to not to use an epsilon adjusted p-value for the Beta approximation.  
  	  	When TRUE, this prevents observed p-values of 0.  The default is TRUE.}
  	  	
  	   \item{scaleXY}{A boolean indicating whether or to not to scale x and y.  The default is TRUE.}
}
\value{
An object with the corresponding GSEA results.  If approx="norm" an \code{npGSEAResultNorm} 
object is returned.  If approx="beta" a \code{npGSEAResultBeta} object is returned.  If approx="chiSq" 
a \code{npGSEAResultChiSq object} is returned.  If set is a \code{GeneSetCollection} (i.e., multiple 
sets of interest), then the corresponding \code{npGSEAResultNormCollection}, 
\code{npGSEAResultBetaCollection}, or \code{npGSEAResultChiSqCollection} is returned.
  }
\description{
  This function calculates the permutation gene set enrichment analysis test statistic and p-value 
  without actually running the permutation.  We account for the covariance among the genes within 
  the set and approximate the corresponding permutation distribution.  For more details on the 
  method see Larson and Owen (2014).
}
\author{
  Jessica L. Larson and Art Owen
}
\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15)
}