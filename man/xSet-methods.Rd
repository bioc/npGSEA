\name{xSet-methods}
\docType{methods}
\alias{xSet}
\alias{xSet-methods}
\alias{xSet,npGSEAResultNorm-method}
\alias{xSet,npGSEAResultBeta-method}
\alias{xSet,npGSEAResultChiSq-method}
\alias{xSet,npGSEAResultNormCollection-method}
\alias{xSet,npGSEAResultBetaCollection-method}
\alias{xSet,npGSEAResultChiSqCollection-method}
\title{ ~~ Methods for Function \code{xSet}  ~~}
\description{
This function returns the scaled and centered expression data 
for all genes in the corresponding \code{GeneSet} in the given experiment, 
or a list of such matrices for each set in a \code{GeneSetCollection}.  
This method is applicable for all three approximation methods.
}

\usage{
    xSet(object)
}
\arguments{
    \item{object}{An object of type \code{npGSEAResultNorm}, \code{npGSEAResultBeta}, 
    \code{npGSEAResultChiSq},\code{npGSEAResultNormCollection}, 
    \code{npGSEAResultBetaCollection}, or \code{npGSEAResultChiSqCollection} }
}

\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Returns the centered and 
    scaled X matrix used for analysis to create a \code{npGSEAResultNorm} object}  
    \item{\code{signature(object = "npGSEAResultBeta")}}{Returns the centered and 
    scaled X matrix used for analysis to create a \code{npGSEAResultBeta} object}    
    \item{\code{signature(object = "npGSEAResultChiSq")}}{Returns the centered and 
    scaled X matrix used for analysis to create a \code{npGSEAResultChiSq} object}
    \item{\code{signature(object = "npGSEAResultNormCollection")}}{Returns a list of 
    centered and scaled X matrices used for analysis to create the 
    \code{npGSEAResultNormCollection} objects (1 for each set)}    
    \item{\code{signature(object = "npGSEAResultBetaCollection")}}{Returns a list of 
    centered and scaled X matrices used for analysis to create the 
    \code{npGSEAResultBetaCollection} objects (1 for each set)}    
    \item{\code{signature(object = "npGSEAResultChiSqCollection")}}{Returns a list of 
    centered and scaled X matrices used for analysis to create the 
    \code{npGSEAResultChiSqCollection} objects (1 for each set)}
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
    head( xSet(res) )
}


\keyword{methods}

