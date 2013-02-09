\name{copynumbR.read.segmented}
\alias{copynumbR.read.segmented}
\title{copynumbR read segmented}
\description{Read segmented data and turn it into an ExpressionSet.}
\usage{copynumbR.read.segmented(filename, clinical, gene = FALSE, mad = 0, 
    geneMap = NULL, genes.as.features = TRUE, ...)}
\arguments{
  \item{filename}{The filename of segmented data.}
  \item{clinical}{A data frame with clinical annotation for the phenoData slot of the output
ExpressionSet.}
  \item{gene}{Either segments or genes as features. }
  \item{mad}{Filter low-variance segments, typically useful for clustering.
Calculates the mean absolute deviation across samples
for each rows and drops rows that are not above the percentile defined.}
  \item{geneMap}{The gene mapping if gene==TRUE. See ?getRS. Maps for hg17, hg18 and hg19
are provided by this package (geneMap_hgXX).}
  \item{genes.as.features}{If gene==TRUE, use genes as feature names. Note that feature names
have to be unique and only the first occurence of a gene in the geneMap
is used if genes.as.features==TRUE.}
  \item{\dots}{Additional arguments passed to the copynumbR.eset function.}
}

\value{An ExpressionSet.    }

\author{Markus Riester}




\examples{
library(copynumbR)

clinical <- read.csv(system.file("extdata", 
    "stransky_bladder_clinical.csv", package="copynumbR"))
eset <- copynumbR.read.segmented(system.file("extdata", 
    "stransky_bladder.glad", package="copynumbR"), clinical)

# Plot the distribution of copy numbers of segments in all samples
plot(exprs(eset))

eset.genes <- copynumbR.read.segmented(system.file("extdata",
 "stransky_bladder.glad", package="copynumbR"), clinical, gene=TRUE,
 geneMap=geneMap_hg17)

boxplot(exprs(eset.genes)["MYC",])
}
