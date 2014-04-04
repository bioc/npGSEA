\name{getIncidence}
\alias{getIncidence}
\title{Determines the incidence of a gene set in a list of genes.}
\usage{
  getIncidence(universeIDs, set)
}
\arguments{
    \item{universeIDs}{A vector containing the list of 
        possible gene ids in the universe (or experiment).}

    \item{set}{A GeneSet object containing a set of genes of
        interest}
}
\value{
    A list of inSet and inExp.  inSet is a vector with the same length as universeIDs.  
    Each value of inSet is 1 if the gene is in the set and 0 otherwise.  
    inExp is a vector with the same length as geneIds(set), the number of genes in the set.  
    Each value of inExp is 1 if the gene is in universeIDs and 0 otherwise.
}
\description{
    Calculates the incidence of a gene set in an experiment
}
\details{
    getIncidence returns an incidence vector of the location of the genes within a gene set in a list of genes in an experiment and vise-versa.
}
\examples{
    geneSetABC15 <- GeneSet(geneIds=letters[1:15], setName="setABC15")
    getIncidence(letters, geneSetABC15)
}
\author{
    Jessica L. Larson and Art Owen
}

