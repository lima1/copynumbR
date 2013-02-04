\name{copynumbR.gistic.lesions}
\alias{copynumbR.gistic.lesions}
\title{copynumbR gistic lesions}
\description{Read GISTIC all_lesions file}
\usage{copynumbR.gistic.lesions(filename, clinical, data.col = 10, ...)}
\arguments{
  \item{filename}{The filename of the GISTIC output file}
  \item{clinical}{A data frame with clinical annotation for the phenoData slot of the output
ExpressionSet}
  \item{data.col}{Start column of the sample data, no need to change unless GISTIC output
changes }
  \item{\dots}{Additional arguments passed to copynumbR.eset()}
}

\value{ExpressionSet containing the significant GISTIC output    }

\author{Markus Riester}





