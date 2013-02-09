\name{copynumbR.read.segmented}
\alias{copynumbR.read.segmented}
\title{copynumbR read segmented}
\description{Read segmented data and turn it into an ExpressionSet}
\usage{copynumbR.read.segmented(filename, clinical, gene = FALSE, mad = 0, 
    geneMap = NULL, ...)}
\arguments{
  \item{filename}{The filename of segmented data}
  \item{clinical}{A data frame with clinical annotation for the phenoData slot of the output
ExpressionSet}
  \item{gene}{Either segments or genes as features. }
  \item{mad}{filter low-variance segments, typically useful for clustering}
  \item{geneMap}{The gene mapping if gene==TRUE. See ?getRS.}
  \item{\dots}{Additional arguments passed to the copynumbR.eset function}
}



\author{Markus Riester}




\examples{
library(copynumbR)
clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
eset <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)
# Plot the distribution of copy numbers of segments in all samples
plot(exprs(eset))

eset.genes <- copynumbR.read.segmented(system.file("extdata",
 "stransky_bladder.glad", package="copynumbR"), clinical, gene=TRUE,
 geneMap=geneMap_hg17)

boxplot(exprs(eset.genes)["MYC",])
}
