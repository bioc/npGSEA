\name{npGSEAResultChiSq-class}
\Rdversion{1.1}
\docType{class}
\alias{npGSEAResultChiSq-class}

\title{Class \code{"npGSEAResultChiSq"}}
\description{
Objects of this class store results from running npGSEA with the Chi-square approximation. 
}
\section{Objects from the Class}{
Objects can be created by calls of \code{npGSEA}.
}
\section{Slots}{
    \describe{
        \item{\code{geneSetName}:}{Object of class \code{"character"}, the name of the geneSet }
        \item{\code{chiSqStat}:}{Object of class \code{"numeric"} , the test statistic, scaled to the corresponding chi-sq distribution }
        \item{\code{ChatGw}:}{Object of class \code{"numeric"}, the test statistic for the set }
        \item{\code{sigmaSq}:}{Object of class \code{"numeric"}, the variance }
        \item{\code{DF}:}{Object of class \code{"numeric"}, the degrees of freedom }
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
showClass("npGSEAResultChiSq")
}
\keyword{classes}
