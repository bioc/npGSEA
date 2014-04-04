\name{DF-methods}
\docType{methods}
\alias{DF}
\alias{DF-methods}
\alias{DF,npGSEAResultChiSq-method}
\alias{DF,npGSEAResultChiSqCollection-method}
\title{ ~~ Methods for Function \code{DF}  ~~}
\description{
This function returns the degrees of freedom for the chi-sq approximation of a 
corresponding \code{GeneSet} 
or a list of degrees of freedom for a \code{GeneSetCollection}.  
This method is applicable for only the chi-sq approximation method.
}
\usage{DF(object)}
\arguments{
    \item{object}{An object of type 
        \code{npGSEAResultChiSq} or \code{npGSEAResultChiSqCollection} }
}
\section{Methods}{
\describe{
        \item{\code{signature(object = "npGSEAResultChiSq")}}{Returns the degrees of 
        freedom for a \code{npGSEAResultChiSq} object}
        \item{\code{signature(object = "npGSEAResultChiSqCollection")}}{Returns a 
        list of the degrees of freedom for a \code{npGSEAResultChiSqCollection} 
        objects (1 for each set)}
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
    DF(res)
}
\author{Jessica L. Larson}
\keyword{methods}
