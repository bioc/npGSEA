\name{pTwoSided-methods}
\docType{methods}
\alias{pTwoSided}
\alias{pTwoSided-methods}
\alias{pTwoSided,npGSEAResultNorm-method}
\alias{pTwoSided,npGSEAResultBeta-method}
\alias{pTwoSided,npGSEAResultChiSq-method}
\alias{pTwoSided,npGSEAResultNormCollection-method}
\alias{pTwoSided,npGSEAResultBetaCollection-method}
\alias{pTwoSided,npGSEAResultChiSqCollection-method}
\title{ ~~ Methods for Function \code{pTwoSided}  ~~}
\description{
This function returns the two-sided p-value for the corresponding \code{GeneSet} 
or a list of p-values for a \code{GeneSetCollection}.  
This method is applicable for all three approximation methods.
}
\usage{pTwoSided(object)}
\arguments{
  \item{object}{An object of type \code{npGSEAResultNorm}, \code{npGSEAResultBeta}, 
    \code{npGSEAResultChiSq}, \code{npGSEAResultNormCollection}, 
    \code{npGSEAResultBetaCollection}, or \code{npGSEAResultChiSqCollection} }
}
\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Returns a two-sided 
    p-value for a \code{npGSEAResultNorm} object}
    \item{\code{signature(object = "npGSEAResultBeta")}}{Returns a two-sided 
    p-value for a \code{npGSEAResultBeta} object}    
    \item{\code{signature(object = "npGSEAResultChiSq")}}{Returns a two-sided 
    p-value for a \code{npGSEAResultChiSq} object}
    \item{\code{signature(object = "npGSEAResultNormCollection")}}{Returns a list 
    of left-sided p-values for a \code{npGSEAResultNormCollection} objects 
    (1 for each set)}
    \item{\code{signature(object = "npGSEAResultBetaCollection")}}{Returns a list of 
    left-sided p-values for a \code{npGSEAResultBetaCollection} objects (1 for each set)} 
    \item{\code{signature(object = "npGSEAResultChiSqCollection")}}{Returns a list of 
    two-sided p-values for a \code{npGSEAResultChiSqCollection} objects (1 for each set)}
    }
}

\seealso{
    \code{\linkS4class{npGSEAResultNorm}}-class, \code{pRight}
}

\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15)
    pTwoSided (res)
}
\author{
Jessica L. Larson
}
\keyword{methods}
