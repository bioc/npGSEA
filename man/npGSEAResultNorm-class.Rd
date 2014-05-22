\name{npGSEAResultNorm-class}
\Rdversion{1.1}
\docType{class}
\alias{npGSEAResultNorm-class}

\title{Class \code{"npGSEAResultNorm"}}
\description{
Objects of this class store results from running npGSEA with the Gaussian approximation. 
}
\section{Objects from the Class}{
Objects can be created by calls of \code{npGSEA}.
}
\section{Slots}{
    \describe{
        \item{\code{geneSetName}:}{Object of class \code{"character"}, the name of the geneSet }
        \item{\code{zStat}:}{Object of class \code{"numeric"} , the test statistic, scaled to a standard normal }
        \item{\code{ThatGw}:}{Object of class \code{"numeric"}, the test statistic for the set }
        \item{\code{varThatGw}:}{Object of class \code{"numeric"}, the variance of ThatGw }
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
\code{\linkS4class{npGSEAResultChiSq}}
}
\examples{
showClass("npGSEAResultNorm")
}
\keyword{classes}
