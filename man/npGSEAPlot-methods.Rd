\name{npGSEAPlot-methods}
\docType{methods}
\alias{npGSEAPlot}
\alias{npGSEAPlot-methods}
\alias{npGSEAPlot,npGSEAResultNorm-method}
\alias{npGSEAPlot,npGSEAResultBeta-method}
\alias{npGSEAPlot,npGSEAResultChiSq-method}
\title{ ~~ Methods for Function \code{npGSEAPlot}  ~~}
\description{
This function plots the reference distribution and the 
corresponding scaled statistic (Z, Beta, or Chi-sq) from the npGSEA analysis 
 for a given \code{GeneSet}.  
This method is applicable for all three approximation methods.
}
\usage{
    npGSEAPlot(object)
}
\arguments{
    \item{object}{ An object of type \code{npGSEAResultNorm}, \code{npGSEAResultBeta}, 
      or \code{npGSEAResultChiSq} }
}
\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Plots the Z-statistic for a 
    \code{npGSEAResultNorm} object and the standard normal distribution}
    \item{\code{signature(object = "npGSEAResultBeta")}}{Plots the beta statistic 
    for a \code{npGSEAResultBeta} object and the corresponding reference beta 
    distribution (with alpha and beta calculated from \code{npGSEA}).}   
    \item{\code{signature(object = "npGSEAResultChiSq")}}{Plots the beta 
    statistic for a \code{npGSEAResultChiSq} object and the corresponding reference 
    chi-squared distribution (with degrees of freedom calculated from \code{npGSEA}).}
}
}
\seealso{\code{\linkS4class{npGSEAResultNorm}}-class}
\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15)
    ##npGSEAPlot (res)
}
\author{
Jessica L. Larson
}
\keyword{methods}
