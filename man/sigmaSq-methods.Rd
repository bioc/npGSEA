\name{sigmaSq-methods}
\docType{methods}
\alias{sigmaSq}
\alias{sigmaSq-methods}
\alias{sigmaSq,npGSEAResultNorm-method}
\alias{sigmaSq,npGSEAResultBeta-method}
\alias{sigmaSq,npGSEAResultChiSq-method}
\alias{sigmaSq,npGSEAResultNormCollection-method}
\alias{sigmaSq,npGSEAResultBetaCollection-method}
\alias{sigmaSq,npGSEAResultChiSqCollection-method}
\title{ ~~ Methods for Function \code{sigmaSq}  ~~}
\description{
This function returns the corresponding variance of the statistic (linear or quadratic) from the npGSEA analysis 
 for a given \code{GeneSet},
or a list of these variances for a given \code{GeneSetCollection}.  
This method is applicable for all three approximation methods.
}
\usage{sigmaSq(object)}
\arguments{
  \item{object}{An object of type \code{npGSEAResultNorm}, \code{npGSEAResultBeta}, 
      \code{npGSEAResultChiSq},
      \code{npGSEAResultNormCollection}, \code{npGSEAResultBetaCollection}, or
       \code{npGSEAResultChiSqCollection} }
}

\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Returns the variance 
    of the linear statistic for a \code{npGSEAResultNorm} object}  
    \item{\code{signature(object = "npGSEAResultBeta")}}{Returns the variance 
    of the linear statistic for a \code{npGSEAResultBeta} object} 
    \item{\code{signature(object = "npGSEAResultChiSq")}}{Returns the variance 
    of the quadratic statistic for a \code{npGSEAResultChiSq} object} 
    \item{\code{signature(object = "npGSEAResultNormCollection")}}{Returns a list 
    of the variances of the linear statistics for a \code{npGSEAResultNormCollection} 
    objects (1 for each set)}
    \item{\code{signature(object = "npGSEAResultBetaCollection")}}{Returns a list of 
    the variances of the linear statistics for a \code{npGSEAResultBetaCollection} 
    objects (1 for each set)}
    \item{\code{signature(object = "npGSEAResultChiSqCollection")}}{Returns a list of 
    the variances of the linear statistics for a \code{npGSEAResultChiSqCollection} 
    objects (1 for each set)} 
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
    sigmaSq(res)
}
\author{
Jessica L. Larson
}
\keyword{methods}
