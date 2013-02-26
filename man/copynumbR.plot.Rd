\name{copynumbR.plot}
\alias{copynumbR.plot}
\title{copynumbR plot}
\description{Chromosome frequency plot.}
\usage{copynumbR.plot(esets, window = 5e+05, labels = paste(names(esets), 
    " (n = ", sapply(esets, ncol), ")", sep = ""), gain = 0.1, 
    loss = -0.1, sma = 10, from.chr = 1, to.chr = 22, ylab = "Loss / Gain", 
    xlab = "Chromosome", font.size = 12, centromere.file = "hg18")}
\arguments{
  \item{esets}{List of ExpressionSets with segmented copy number data.}
  \item{window}{Window size in base pairs. Copy numbers will be the median in each window.}
  \item{labels}{The labels of the ExpressionSets.}
  \item{gain}{Minimum log2 ratio for a copy number gain.}
  \item{loss}{Maximum log2 ratio for a copy number loss.}
  \item{sma}{Smooth copy numbers with a simple moving average.}
  \item{from.chr}{Start plotting at this chromosome.}
  \item{to.chr}{End plotting at this chromosome.}
  \item{ylab}{The y-axis label.}
  \item{xlab}{The x-axis label.}
  \item{font.size}{The font size.}
  \item{centromere.file}{File containing the centromere locations for each chromosome.
These files are already provided for hg17-hg19.
Select these with centromere.file=c("hg17", "hg18","hg19").}
}

\value{A ggplot2 object.}

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

# find the non-muscle-invasive samples
idx.noninvasive <- grepl("Ta|T1", eset$T)

# create a list of ExpressionSets with the low- and high stage samples
eset.stage <- list("Non-Muscle-Invasive"=eset[,idx.noninvasive],
    "Invasive"=eset[,!idx.noninvasive])

# now compare the copy numbers of these two groups
p <- copynumbR.plot(eset.stage, centromere.file="hg17", sma=0)

plot(p) 
}
