\name{copynumbR.plot}
\alias{copynumbR.plot}
\title{copynumbR plot}
\description{Chromosome frequency plot}
\usage{copynumbR.plot(esets, window = 5e+05, labels = paste(names(esets), 
    " (n = ", sapply(esets, ncol), ")", sep = ""), gain = 0.1, 
    loss = -0.1, sma = 10, from.chr = 1, to.chr = 22, ylab = "Loss (\%) / Gain (\%)", 
    xlab = "Chromosome", centromer.file = system.file("extdata", 
        "hg18centromer.txt", package = "copynumbR"))}
\arguments{
  \item{esets}{List of ExpressionSets}
  \item{window}{Window size in base pairs. Copy numbers will be averaged in each window}
  \item{labels}{The labels of the ExpressionSets.}
  \item{gain}{Minimum log2 ratio for a copy number gain}
  \item{loss}{Maximum log2 ratio for a copy number loss}
  \item{sma}{Smooth copy numbers with a simple moving average}
  \item{from.chr}{Start plotting at this chromosome}
  \item{to.chr}{End plotting at this chromosome}
  \item{ylab}{The y-axis label}
  \item{xlab}{The x-axis label}
  \item{centromer.file}{
}
}

\value{A ggplot2 object}

\author{Markus Riester}




\examples{
library(copynumbR)
clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
eset <- copynumbR.read.segmented(system.file("extdata", "stransky_bladder.glad", package="copynumbR"), clinical)
p <- copynumbR.plot(list(eset), labels="Stransky")
plot(p) 
}
