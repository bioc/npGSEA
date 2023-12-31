\name{npGSEAResultChiSqCollection-class}
\Rdversion{1.1}
\docType{class}
\alias{npGSEAResultChiSqCollection-class}

\title{Class \code{"npGSEAResultChiSqCollection"}}
\description{
Objects of this class store results from running npGSEA with the Chi-square approximation with a GeneSetCollection.   
npGSEAResultChiSqCollection objects contain a list of npGSEAResultChiSq objects (one result for each GeneSet).
}
\section{Objects from the Class}{
Objects can be created by calls of \code{npGSEA}.
}
\author{
Jessica L. Larson
}

\seealso{
\code{\linkS4class{npGSEAResultNormCollection}}
}
\examples{
showClass("npGSEAResultChiSqCollection")
}
\keyword{classes}
\references{
Jessica L Larson and Art B Owen: Moment based gene set tests. BMC Bioinformatics 2015, 16:132. \url{http://www.biomedcentral.com/1471-2105/16/132}
}
