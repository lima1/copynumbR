\name{copynumbR.gistic.genes.band}
\alias{copynumbR.gistic.genes.band}
\title{copynumbR gistic genes band}
\description{Extract the GISTIC target genes }
\usage{copynumbR.gistic.genes.band(filename, band)}
\arguments{
  \item{filename}{The GISTIC output file amp_genes.conf_95.txt}
  \item{band}{The chromosome band}
}

\value{Genes in the specified chromosome band    }

\author{Markus Riester}




\examples{
library(copynumbR)
# extract the gene symbols in the GISTIC peak 8q24.21
band <-
copynumbR.gistic.genes.band(system.file("extdata/gistic_stransky_bladder",
    "amp_genes.conf_95.txt", package="copynumbR"), "8q24.21")
}
