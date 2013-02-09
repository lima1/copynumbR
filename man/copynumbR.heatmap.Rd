\name{copynumbR.heatmap}
\alias{copynumbR.heatmap}
\title{copynumbR heatmap}
\description{Chromosome Heatmap}
\usage{copynumbR.heatmap(eset, window = 1e+06, gain = 0.1, loss = -0.1, 
    Colv = NULL, distfun = spearman.dist, hclustfun = hclust, 
    from.chr = 1, to.chr = 22, start = NULL, end = NULL, ensembl = NULL, 
    plot = TRUE, hide.labels = FALSE, ylab = ("Chromosome"), 
    centromere.file = "hg18")}
\arguments{
  \item{eset}{ExpressionSet, typically created with copynumbR.read.segmented}
  \item{window}{Window size in base pairs. Copy numbers will the median in each window}
  \item{gain}{Minimum log2 ratio for a copy number gain}
  \item{loss}{Maximum log2 ratio for a copy number loss}
  \item{Colv}{Cluster columns.}
  \item{distfun}{The distance function if Colv is NULL}
  \item{hclustfun}{The cluster function if Colv is NULL}
  \item{from.chr}{Start plotting at this chromosome}
  \item{to.chr}{End plotting at this chromosome}
  \item{start}{
}
  \item{end}{
}
  \item{ensembl}{
}
  \item{plot}{Plot the heatmap, otherwise the dendrogram and heatmap are returned as
list of ggplot2 objects}
  \item{hide.labels}{Hide the sample labels on the bottom}
  \item{ylab}{The y-axis label}
  \item{centromere.file}{File containing the centromere locations for each chromosome.
These files are already provided for hg17-hg19
curl -s
"http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz" |
gunzip -c | grep acen
Select these with centromere.file=c("hg17", "hg18","hg19")}
}

\value{List of two ggplot2 objects (heatmap and dendrogram)}

\author{Markus Riester}




\examples{
library(copynumbR)

clinical <- read.csv(system.file("extdata", 
    "stransky_bladder_clinical.csv", package="copynumbR"))

eset <- copynumbR.read.segmented(system.file("extdata", 
    "stransky_bladder.glad", package="copynumbR"), clinical)

p <- copynumbR.heatmap(eset, centromere.file="hg17")
}
