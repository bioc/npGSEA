\name{betaHats-methods}
\docType{methods}
\alias{betaHats}
\alias{betaHats-methods}
\alias{betaHats,npGSEAResultNorm-method}
\alias{betaHats,npGSEAResultBeta-method}
\alias{betaHats,npGSEAResultChiSq-method}
\alias{betaHats,npGSEAResultNormCollection-method}
\alias{betaHats,npGSEAResultBetaCollection-method}
\alias{betaHats,npGSEAResultChiSqCollection-method}
\title{ ~~ Methods for Function \code{betaHats}  ~~}
\description{
This function returns the betaHats 
for all genes in the corresponding \code{GeneSet} in the given experiment, 
or a list of such vectors for each set in a \code{GeneSetCollection}.  
This corresponds to each gene's contrubution to the test statistic.
This method is applicable for all three approximation methods.
}

\usage{
    betaHats(object)
}
\arguments{
    \item{object}{An object of type \code{npGSEAResultNorm}, \code{npGSEAResultBeta}, 
    \code{npGSEAResultChiSq},\code{npGSEAResultNormCollection}, 
    \code{npGSEAResultBetaCollection}, or \code{npGSEAResultChiSqCollection} }
}

\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Returns the betaHats 
    	used for analysis to create a \code{npGSEAResultNorm} object}  
    \item{\code{signature(object = "npGSEAResultBeta")}}{Returns the betaHats 
    used for analysis to create a \code{npGSEAResultBeta} object}    
    \item{\code{signature(object = "npGSEAResultChiSq")}}{Returns the betaHats
    used for analysis to create a \code{npGSEAResultChiSq} object}
    \item{\code{signature(object = "npGSEAResultNormCollection")}}{Returns the betaHats 
    	used for analysis to create the 
    \code{npGSEAResultNormCollection} objects (1 for each set)}    
    \item{\code{signature(object = "npGSEAResultBetaCollection")}}{Returns a list of 
    the betaHats used for analysis to create the 
    \code{npGSEAResultBetaCollection} objects (1 for each set)}    
    \item{\code{signature(object = "npGSEAResultChiSqCollection")}}{Returns a list of 
    the betaHats used for analysis to create the 
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
    betaHats(res) 
}


\keyword{methods}

