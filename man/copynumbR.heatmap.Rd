\name{copynumbR.heatmap}
\alias{copynumbR.heatmap}
\title{copynumbR heatmap}
\description{Plot a chromosome heatmap.}
\usage{copynumbR.heatmap(eset, window = 1e+06, gain = 0.1, loss = -0.1, 
    Colv = NULL, distfun = spearman.dist, hclustfun = hclust, 
    from.chr = 1, to.chr = 22, start = NULL, end = NULL, mart = NULL, 
    plot = TRUE, hide.labels = FALSE, ylab = ("Chromosome"), 
    centromere.file = "hg18")}
\arguments{
  \item{eset}{An ExpressionSet, typically created with copynumbR.read.segmented().}
  \item{window}{Window size in base pairs. Copy numbers will be the median in each window.}
  \item{gain}{Minimum log2 ratio for a copy number gain.}
  \item{loss}{Maximum log2 ratio for a copy number loss.}
  \item{Colv}{Hierarchical clustering of class hclust of the columns. 
If NULL, then calculated as specified by distfun and hclustfun.
If NA, no ordering of columns is done.}
  \item{distfun}{The distance function if Colv is NULL.}
  \item{hclustfun}{The cluster function if Colv is NULL.}
  \item{from.chr}{Start plotting at this chromosome.}
  \item{to.chr}{End plotting at this chromosome.}
  \item{start}{Zoom in, start from the specified coordinates (bp) of from.chr}
  \item{end}{Zoom in, stop at the specified coordinates (bp) of to.chr. Only
useful when from.chr and to.chr are equal}
  \item{mart}{A biomaRt object, used for a zoom in to annotate genes.}
  \item{plot}{Plot the heatmap, otherwise the dendrogram and heatmap are returned as
list of ggplot2 objects.}
  \item{hide.labels}{Hide the sample labels on the bottom.}
  \item{ylab}{The y-axis label.}
  \item{centromere.file}{File containing the centromere locations for each chromosome.
These files are already provided for hg17-hg19.
Select these with centromere.file=c("hg17", "hg18","hg19").}
}

\value{List of two ggplot2 objects (heatmap and dendrogram).}

\author{Markus Riester}
\note{The centromere files were generated with: \code{curl -s
"http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz" |
gunzip -c | grep acen}}



\examples{
library(copynumbR)

clinical <- read.csv(system.file("extdata", 
    "stransky_bladder_clinical.csv", package="copynumbR"))

eset <- copynumbR.read.segmented(system.file("extdata", 
    "stransky_bladder.glad", package="copynumbR"), clinical)

p <- copynumbR.heatmap(eset, centromere.file="hg17")
}
