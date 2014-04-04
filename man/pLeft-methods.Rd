\name{pLeft-methods}
\docType{methods}
\alias{pLeft}
\alias{pLeft-methods}
\alias{pLeft,npGSEAResultNorm-method}
\alias{pLeft,npGSEAResultBeta-method}
\alias{pLeft,npGSEAResultNormCollection-method}
\alias{pLeft,npGSEAResultBetaCollection-method}
\title{ ~~ Methods for Function \code{pLeft}  ~~}
\description{
This function returns the left-sided p-value for the corresponding \code{GeneSet} 
or a list of p-values for a \code{GeneSetCollection}.  
This method is only applicable for the normal and beta approximation methods.
}
\usage{
pLeft(object)
}
\arguments{
\item{object}{An object of type \code{npGSEAResultNorm}, \code{npGSEAResultBeta},
    \code{npGSEAResultNormCollection}, or \code{npGSEAResultBetaCollection} }
}
\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Returns a left-sided p-value 
    for a \code{npGSEAResultNorm} object}
    \item{\code{signature(object = "npGSEAResultBeta")}}{Returns a left-sided p-value 
    for a \code{npGSEAResultBeta} object}
    \item{\code{signature(object = "npGSEAResultNormCollection")}}{Returns a list of 
    left-sided p-values for a \code{npGSEAResultNormCollection} objects (1 for each set)}
    \item{\code{signature(object = "npGSEAResultBetaCollection")}}{Returns a list of 
    left-sided p-values for a \code{npGSEAResultBetaCollection} objects (1 for each set)}
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
    pLeft (res)
}
\author{
Jessica L. Larson
}
\keyword{methods}
