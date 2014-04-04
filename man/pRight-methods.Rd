\name{pRight-methods}
\docType{methods}
\alias{pRight}
\alias{pRight-methods}
\alias{pRight,npGSEAResultNorm-method}
\alias{pRight,npGSEAResultBeta-method}
\alias{pRight,npGSEAResultNormCollection-method}
\alias{pRight,npGSEAResultBetaCollection-method}
\title{ ~~ Methods for Function \code{pRight}  ~~}
\description{
This function returns the right-sided p-value for the corresponding \code{GeneSet} 
or a list of p-values for a \code{GeneSetCollection}.  
This method is only applicable for the normal and beta approximation methods.
}
\usage{pRight(object)}
\arguments{
\item{object}{An object of type \code{npGSEAResultNorm}, \code{npGSEAResultBeta}, 
\code{npGSEAResultNormCollection}, or \code{npGSEAResultBetaCollection} }
}
\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Returns a right-sided 
    p-value for a \code{npGSEAResultNorm} object}
    \item{\code{signature(object = "npGSEAResultBeta")}}{Returns a right-sided 
    p-value for a \code{npGSEAResultBeta} object}
    \item{\code{signature(object = "npGSEAResultNormCollection")}}{Returns a list 
    of right-sided p-values for a \code{npGSEAResultNormCollection} objects 
    (1 for each set)}
    \item{\code{signature(object = "npGSEAResultBetaCollection")}}{Returns a list of 
    right-sided p-values for a \code{npGSEAResultBetaCollection} objects 
    (1 for each set)}
}
}
\seealso{
\code{\linkS4class{npGSEAResultNorm}}-class, \code{pLeft}
}
\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15)
    pRight (res)
}
\author{
Jessica L. Larson
}
\keyword{methods}
