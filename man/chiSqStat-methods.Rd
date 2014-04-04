\name{chiSqStat-methods}
\docType{methods}
\alias{chiSqStat}
\alias{chiSqStat-methods}
\alias{chiSqStat,npGSEAResultChiSq-method}
\alias{chiSqStat,npGSEAResultChiSqCollection-method}
\title{ ~~ Methods for Function \code{chiSqStat}  ~~}
\description{
This function returns the chi-sq statistic (which is compared to a reference 
Chi-sq distribution) for the chi-sq approximation of \code{npGSEA} for 
a corresponding \code{GeneSet} 
or a list of these statistics for a \code{GeneSetCollection}.  
This method is applicable for only the chi-sq approximation method.
}
\usage{chiSqStat(object)}
\arguments{
    \item{object}{An object of type 
        \code{npGSEAResultChiSq} or \code{npGSEAResultChiSqCollection} }
}
\section{Methods}{
\describe{
        \item{\code{signature(object = "npGSEAResultChiSq")}}{Returns the chi-sq statistic for a 
        \code{npGSEAResultChiSq} object}
        \item{\code{signature(object = "npGSEAResultChiSqCollection")}}{Returns a list of 
        the chi-sq statistics for a \code{npGSEAResultChiSqCollection} objects (1 for each set)}
        }
}
\seealso{\code{\linkS4class{npGSEAResultChiSq}}-class}
\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15, approx= "chiSq")  
    chiSqStat(res)
}
\author{Jessica L. Larson}
\keyword{methods}
