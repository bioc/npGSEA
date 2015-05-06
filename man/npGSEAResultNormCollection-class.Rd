\name{npGSEAResultNormCollection-class}
\Rdversion{1.1}
\docType{class}
\alias{npGSEAResultNormCollection-class}

\title{Class \code{"npGSEAResultNormCollection"}}
\description{
Objects of this class store results from running npGSEA with the Gaussian approximation with a GeneSetCollection.   
npGSEAResultNormCollection objects contain a list of npGSEAResultNorm objects (one result for each GeneSet).
}
\section{Objects from the Class}{
Objects can be created by calls of \code{npGSEA}.
}
\author{
Jessica L. Larson
}

\seealso{
\code{\linkS4class{npGSEAResultChiSqCollection}}
}
\examples{
showClass("npGSEAResultNormCollection")
}
\keyword{classes}
\references{
Jessica L Larson and Art B Owen: Moment based gene set tests. BMC Bioinformatics 2015, 16:132. \url{http://www.biomedcentral.com/1471-2105/16/132}
}