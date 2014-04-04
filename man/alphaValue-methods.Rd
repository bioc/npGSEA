\name{alphaValue-methods}
\docType{methods}
\alias{alphaValue}
\alias{alphaValue-methods}
\alias{alphaValue,npGSEAResultBeta-method}
\alias{alphaValue,npGSEAResultBetaCollection-method}
\title{ ~~ Methods for Function \code{alphaValue}  ~~}
\description{
This function returns the corresponding alpha value for the reference 
beta distribution for the \code{npGSEA} analysis in the gene set in the 
given experiment.  This method is applicable for only the beta approximation method.
}
\usage{alphaValue(object)}
\arguments{
    \item{object}{An object of type 
        \code{npGSEAResultBeta} or \code{npGSEAResultBetaCollection}}
}
\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultBeta")}}{Returns the value for 
    	alpha for a \code{npGSEAResultBeta} object}
    \item{\code{signature(object = "npGSEAResultBetaCollection")}}{Returns a 
    	list of the alpha values for a \code{npGSEAResultBetaCollection} 
        objects (1 for each set)}}
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
    alphaValue(res)
}
\author{Jessica L. Larson}
\keyword{methods}