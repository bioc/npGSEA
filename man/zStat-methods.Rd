\name{zStat-methods}
\docType{methods}
\alias{zStat}
\alias{zStat-methods}
\alias{zStat,npGSEAResultNorm-method}
\alias{zStat,npGSEAResultNormCollection-method}
\title{ ~~ Methods for Function \code{zStat}  ~~}
\description{
This function returns the Z-statistic (which is compared to a reference 
standard normal distribution) for the normal approximation of \code{npGSEA} for 
a corresponding \code{GeneSet} 
or a list of these statistics for a \code{GeneSetCollection}.  
This method is applicable for only the normal approximation method.
}
\usage{
    zStat(object)
}
\arguments{
    \item{object}{An object of type 
        \code{npGSEAResultNorm} or \code{npGSEAResultNormCollection}}
}

\section{Methods}{
\describe{
    \item{\code{signature(object = "npGSEAResultNorm")}}{Returns the 
    Z-statistic for a \code{npGSEAResultNorm} object}
    \item{\code{signature(object = "npGSEAResultNormCollection")}}{Returns a 
    list of the Z- statistics for a \code{npGSEAResultNormCollection} objects (1 for each set)}
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
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15, approx= "norm")  
    zStat(res)
}
\author{Jessica L. Larson}
\keyword{methods}