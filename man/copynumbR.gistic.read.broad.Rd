\name{copynumbR.gistic.read.broad}
\alias{copynumbR.gistic.read.broad}
\title{copynumbR gistic read broad}
\description{Read GISTIC broad_values_by_arm.txt output file.}
\usage{copynumbR.gistic.read.broad(filename = "broad_values_by_arm.txt", 
    clinical, data.col = 2, ...)}
\arguments{
  \item{filename}{The filename of the GISTIC broad_values_by_arm.txt output file.}
  \item{clinical}{A data frame with clinical annotation for the phenoData slot of the output
ExpressionSet.}
  \item{data.col}{Start column of the sample data, no need to change unless GISTIC output
changes.}
  \item{\dots}{Additional arguments passed to copynumbR.eset().}
}



\author{Markus Riester}




\examples{
library(copynumbR)
clinical <- read.csv(system.file("extdata", 
    "stransky_bladder_clinical.csv", package="copynumbR"))
eset <-
copynumbR.gistic.read.broad(system.file("extdata/gistic_stransky_bladder",
    "broad_values_by_arm.txt", package="copynumbR"), clinical)

# display the distributions of copy numbers of all chromosome arms
boxplot(t(exprs(eset)))
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      