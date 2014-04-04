\name{geneSetName-methods}
\docType{methods}
\alias{geneSetName}
\alias{geneSetName-methods}
\alias{geneSetName,npGSEAResultNorm-method}
\alias{geneSetName,npGSEAResultBeta-method}
\alias{geneSetName,npGSEAResultChiSq-method}
\alias{geneSetName,npGSEAResultNormCollection-method}
\alias{geneSetName,npGSEAResultBetaCollection-method}
\alias{geneSetName,npGSEAResultChiSqCollection-method}
\title{ ~~ Methods for Function \code{geneSetName}  ~~}
\description{
This function returns the name of the corresponding \code{GeneSet} 
or a list of names for a \code{GeneSetCollection}.  
This method is applicable for all three approximation methods.
}
\usage{geneSetName(object)}
\arguments{
  \item{object}{An object of type \code{npGSEAResultNorm}, \code{npGSEAResultBeta}, 
      \code{npGSEAResultChiSq},
      \code{npGSEAResultNormCollection}, \code{npGSEAResultBetaCollection}, or
       \code{npGSEAResultChiSqCollection} }
}
\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Returns a the name of 
    the gene set from a \code{npGSEAResultNorm} object}
    \item{\code{signature(object = "npGSEAResultBeta")}}{Returns a the name of the 
    gene set from a \code{npGSEAResultBeta} object}
    \item{\code{signature(object = "npGSEAResultChiSq")}}{Returns a the name of 
    the gene set from a \code{npGSEAResultChiSq} object} 
    \item{\code{signature(object = "npGSEAResultNormCollection")}}{Returns a list of 
    the names of the gene sets from a \code{npGSEAResultNormCollection} objects (1 for each set)} 
    \item{\code{signature(object = "npGSEAResultBetaCollection")}}{Returns a list of 
    the names of the gene sets from a \code{npGSEAResultBetaCollection} objects (1 for each set)}
   \item{\code{signature(object = "npGSEAResultChiSqCollection")}}{Returns a list of the 
   	names of the gene sets from a \code{npGSEAResultChiSqCollection} objects (1 for each set)}	
    }
}
\seealso{
  \code{\linkS4class{npGSEAResultNorm}}-class
}
\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15)
    geneSetName (res)
}
\author{
Jessica L. Larson
}
\keyword{methods}