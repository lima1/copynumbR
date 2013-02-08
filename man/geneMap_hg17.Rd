\name{geneMap_hg17}
\alias{geneMap_hg17}
\docType{data}
\title{
UCSC known genes (hg17)
}
\description{
Provides loci annotations for CNTools getRS or
copynumbR.read.segemented.  Data was downloaded from the UCSC genome browser,
sorted with the bedtools, loci without Entrez id removed and only the first
annotation per gene was used.
}
\usage{
data(geneMap_hg17)
}
\format{
A data.frame with columns chromosome (as number), start, end, Entrez id and gene symbol.
}
\source{
Extracted from the UCSC table browser on Feb 7, 2013.
}
\examples{
data(geneMap_hg17)
head(geneMap_hg17)
}
\keyword{datasets}
