\name{copynumbR.gistic.broad}
\alias{copynumbR.gistic.broad}
\title{copynumbR gistic broad}
\description{Read GISTIC broad_values_by_arm.txt output file}
\usage{copynumbR.gistic.broad(filename = "broad_values_by_arm.txt", 
    clinical, data.col = 2, ...)}
\arguments{
  \item{filename}{The filename of the GISTIC output file}
  \item{clinical}{A data frame with clinical annotation for the phenoData slot of the output
ExpressionSet}
  \item{data.col}{Start column of the sample data, no need to change unless GISTIC output
changes }
  \item{\dots}{Additional arguments passed to copynumbR.eset()}
}



\author{Markus Riester}





