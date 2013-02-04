\name{copynumbR.gistic.genes}
\alias{copynumbR.gistic.genes}
\title{copynumbR gistic genes}
\description{Read GISTIC all_data_by_genes file}
\usage{copynumbR.gistic.genes(filename, clinical, data.col = 4, ...)}
\arguments{
  \item{filename}{The filename of the GISTIC output file}
  \item{clinical}{A data frame with clinical annotation for the phenoData slot of the output
ExpressionSet}
  \item{data.col}{Start column of the sample data, no need to change unless GISTIC output
changes }
  \item{\dots}{Additional arguments passed to copynumbR.eset()}
}

\value{ExpressionSet containing the copy numbers for all genes as estimated by
GISTIC}

\author{Markus Riester}





