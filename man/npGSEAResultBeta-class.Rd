\name{npGSEAResultBeta-class}
\Rdversion{1.1}
\docType{class}
\alias{npGSEAResultBeta-class}

\title{Class \code{"npGSEAResultBeta"}}
\description{
Objects of this class store results from running npGSEA with the beta approximation. 
}
\section{Objects from the Class}{
Objects can be created by calls of \code{npGSEA}.
}
\section{Slots}{
    \describe{
        \item{\code{geneSetName}:}{Object of class \code{"character"}, the name of the geneSet }
        \item{\code{betaStat}:}{Object of class \code{"numeric"} , the test statistic, scaled to the corresponding beta distribution }
        \item{\code{ThatGw}:}{Object of class \code{"numeric"}, the test statistic for the set }
        \item{\code{varThatGw}:}{Object of class \code{"numeric"}, the variance of ThatGw }
        \item{\code{alpha}:}{Object of class \code{"numeric"}, the alpha value }
        \item{\code{beta}:}{Object of class \code{"numeric"}, the beta value}
        \item{\code{pLeft}:}{Object of class \code{"numeric"}, the p-value for the left-side hypothesis }
        \item{\code{pRight}:}{Object of class \code{"numeric"}, the p-value for the right-side hypothesis }
        \item{\code{pTwoSided}:}{Object of class \code{"numeric"}, the p-value for the two-sided hypothesis }
        \item{\code{xSet}:}{Object of class \code{"matrix"}, the centered and scaled x data for this set}
        \item{\code{betaHats}:}{Object of class \code{"vector"}, the betaHats for each gene in this set}
    }
}
\author{
Jessica L. Larson
}

\seealso{
\code{\linkS4class{npGSEAResultNorm}}
}
\examples{
showClass("npGSEAResultBeta")
}
\keyword{classes}
