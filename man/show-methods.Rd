\name{summary-methods}
\docType{methods}
\alias{show}
\alias{show-methods}
\alias{show,npGSEAResultNorm-method}
\alias{show,npGSEAResultBeta-method}
\alias{show,npGSEAResultChiSq-method}
\title{ ~~ Methods for Function \code{show} in Package \pkg{base} ~~}
\description{
 ~~ Methods for function \code{show} in package \pkg{base} ~~
}
\section{Methods}{
\describe{
    \item{
        \code{signature(x = "npGSEAResultNorm")},\code{signature(x = "npGSEAResultBeta")},
        \code{signature(x = "npGSEAResultChiSq")}}{These methods display the corresponding 
        statistics (linear or quadratic)for the npGSEA analysis in the gene set in the given experiment.}
    }
}
\examples{
    set.seed(15)
    yFactor <- as.factor( c(rep("treated", 5), rep("control", 5)) )
    xData <- matrix(data=rnorm(length(letters)*10) ,nrow=length(letters), ncol=10)
    rownames(xData) <- letters
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    res <- npGSEA(x = xData, y = yFactor, set = geneSetABC15)
    res
}
\author{Jessica L. Larson}
\keyword{methods}
