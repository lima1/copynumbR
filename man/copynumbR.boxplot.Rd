\name{copynumbR.boxplot}
\alias{copynumbR.boxplot}
\title{copynumbR boxplot}
\description{A boxplot showing correlation of copy number and expression for matched
data.}
\usage{copynumbR.boxplot(eset.cn, eset.expr, cutoffs = c(-Inf, -1.3, 
    -0.1, 0.1, 0.9, Inf), cutoff.labels = c("Homozyg. Deletion", 
    "Heterozyg. Deletion", "Normal", "Gain", "Amplification"), 
    probesets = NULL, min.samples = 3, sqrt = FALSE, xlab = "Copy Number", 
    ylab = "Expression", outlier.shape = NA, plot = TRUE)}
\arguments{
  \item{eset.cn}{ExpressionSet with copy number data. The featureNames and sampleNames
must correspond to the eset.expr provided next. Typically created with
copynumbR.read.segmented(...,gene=TRUE) or
copynumbR.gistic.read.genes(...).}
  \item{eset.expr}{ExpressionSet with expression data. See above, the featureNames must
correspond to the ones of eset.cn. So for Affymetrix data for example,
probe sets need to be collapsed (for example with the WGCNA package). }
  \item{cutoffs}{Copy number cutoffs.}
  \item{cutoff.labels}{The labels of these cutoffs.}
  \item{probesets}{Display only these genes. If null, show all genes. That default is 
obviously only useful for already filtered ExpressionSets.}
  \item{min.samples}{Minimum number of samples in each cutoff category.}
  \item{sqrt}{Square root transform the data for read counts?}
  \item{xlab}{The label of the x-axis.}
  \item{ylab}{The label of the y-axis.}
  \item{outlier.shape}{Display outliers? Passed to geom_boxplot().}
  \item{plot}{Generate a ggplot2 object? If FALSE, then just return a list of means for
each group and gene.}
}

\value{A list of a ggplot2 object (plot) and the mean values in each group for each gene
as data.frame (means)}

\author{Markus Riester}




\examples{
library(copynumbR)
library(ggplot2)
clinical <- read.csv(system.file("extdata", 
    "stransky_bladder_clinical.csv", package="copynumbR"))

eset.genes <- copynumbR.read.segmented(system.file("extdata",
    "stransky_bladder.glad", package="copynumbR"), clinical, 
    gene=TRUE, geneMap=geneMap_hg17)

# load matched expression data
data(PMID17099711.GPL91_eset)

# get the samples with matching expression data
isc <- intersect(sampleNames(eset.genes),
    sampleNames(PMID17099711.GPL91_eset))

p <- copynumbR.boxplot(eset.genes[,isc], PMID17099711.GPL91_eset[,isc],
    probeset=c("MYC", "ADCY8"))

# Highlight samples
 plot(p$plot+geom_jitter(aes(shape=eset.genes[,id]$GENDER),size=4)+
    scale_shape_discrete(name="Gender"))
}
