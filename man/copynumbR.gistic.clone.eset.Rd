\name{copynumbR.gistic.clone.eset}
\alias{copynumbR.gistic.clone.eset}
\title{copynumbR gistic clone eset}
\description{Extract GISTIC peak copy numbers from segmented data.}
\usage{copynumbR.gistic.clone.eset(eset.gistic, eset.segmented)}
\arguments{
  \item{eset.gistic}{An ExpressionSet with GISTIC peaks read with
copynumbR.gistic.read.lesions().}
  \item{eset.segmented}{An ExpressionSet with segmented data, typically from another cohort,
read with copynumbR.read.segmented. This is useful to compare copy numbers at
GISTIC peaks across cohorts.}
}

\value{An ExpressionSet containing copy numbers from eset.segmented of the GISTIC peaks    
in eset.gistic.}

\author{Markus Riester}




\examples{
library(copynumbR)

clinical <- read.csv(system.file("extdata", 
    "stransky_bladder_clinical.csv", package="copynumbR"))

# Read GISTIC peaks
eset.gistic <-
copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
    "all_lesions.conf_95.txt", package="copynumbR"), clinical)

# Read segmented data, typically from another cohort
eset.segmented <- copynumbR.read.segmented(system.file("extdata", 
    "stransky_bladder.glad", package="copynumbR"), clinical)

# We do not have example data of a second cohort, so we just extract the
# GISTIC peaks from the segmented data
eset.cloned <- copynumbR.gistic.clone.eset(eset.gistic, eset.segmented)

plot(exprs(eset.gistic)[2,],
     exprs(eset.cloned)[2,], xlab="GISTIC", ylab="Segmented",
     main="MYC Locus")
}
