\name{geneMap_hg19}
\alias{geneMap_hg19}
\docType{data}
\title{
UCSC known genes (hg19)
}
\description{
Provides loci annotations for CNTools getRS or
copynumbR.read.segemented.  Data was downloaded from the UCSC genome browser,
sorted with the bedtools, loci without Entrez id removed and only the first
annotation per gene was used.
}
\usage{
data(geneMap_hg19)
}
\format{
A data.frame with columns chromosome (as number), start, end, Entrez id and gene symbol.
}
\source{
Extracted from the UCSC table browser on Feb 7, 2013.
}
\examples{
data(geneMap_hg19)
head(geneMap_hg19)
}
\keyword{datasets}
