\name{copynumbR.glad}
\alias{copynumbR.glad}
\title{copynumbR glad}
\description{Read segmented data and turn it into an ExpressionSet}
\usage{copynumbR.glad(filename, clinical, gene = FALSE, mad = 0.1, ...)}
\arguments{
  \item{filename}{The filename of segmented data}
  \item{clinical}{A data frame with clinical annotation for the phenoData slot of the output
ExpressionSet}
  \item{gene}{Either segments or genes as features. The latter currently only support
hg18}
  \item{mad}{filter low-variance segments, typically useful for clustering}
  \item{\dots}{
}
}



\author{Markus Riester}





