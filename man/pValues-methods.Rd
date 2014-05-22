\name{pValues-methods}
\docType{methods}
\alias{pValues}
\alias{pValues-methods}
\alias{pValues,npGSEAResultNorm-method}
\alias{pValues,npGSEAResultBeta-method}
\alias{pValues,npGSEAResultChiSq-method}
\title{ ~~ Methods for Function \code{pValues} ~~}
\description{
 ~~ Methods for function \code{pValues} ~~
}
\section{Methods}{
\describe{
    \item{
        \code{signature(x = "npGSEAResultNorm")},\code{signature(x = "npGSEAResultBeta")},
        \code{signature(x = "npGSEAResultChiSq")}}{These methods display the corresponding 
        p-values for the npGSEA analysis in the gene set in the given experiment.}
    }
}
\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15)
    pValues(res)
}
\author{Jessica L. Larson}
\keyword{methods}
