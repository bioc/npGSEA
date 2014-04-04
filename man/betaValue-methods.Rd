\name{betaValue-methods}
\docType{methods}
\alias{betaValue}
\alias{betaValue-methods}
\alias{betaValue,npGSEAResultBeta-method}
\alias{betaValue,npGSEAResultBetaCollection-method}
\title{ ~~ Methods for Function \code{betaValue}  ~~}
\description{
This function returns the corresponding beta value for the reference 
beta distribution for the \code{npGSEA} analysis in the gene set in the 
given experiment.  This method is applicable for only the beta approximation method.
}
\usage{betaValue(object)}
\arguments{
    \item{object}{An object of type 
        \code{npGSEAResultBeta} or \code{npGSEAResultBetaCollection}}
}
\section{Methods}{
\describe{
        \item{\code{signature(object = "npGSEAResultBeta")} }{Returns the 
        value for beta for a \code{npGSEAResultBeta} object}
        \item{\code{signature(object = "npGSEAResultBetaCollection")} }{Returns a 
        list of the beta values for a \code{npGSEAResultBetaCollection} objects (1 for each set)}
        }
}
\seealso{
    \code{\linkS4class{npGSEAResultBeta}}-class
}
\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15, approx= "beta")  
    betaValue(res)
}
\author{Jessica L. Larson}
\keyword{methods}