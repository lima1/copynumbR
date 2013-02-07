\name{copynumbR.gistic.read.lesions}
\alias{copynumbR.gistic.read.lesions}
\title{copynumbR gistic read lesions}
\description{Read GISTIC all_lesions file}
\usage{copynumbR.gistic.read.lesions(filename, clinical, data.col = 10, 
    ...)}
\arguments{
  \item{filename}{The filename of the GISTIC output file}
  \item{clinical}{A data frame with clinical annotation for the phenoData slot of the output
ExpressionSet}
  \item{data.col}{Start column of the sample data, no need to change unless GISTIC output
changes }
  \item{\dots}{Additional arguments passed to copynumbR.eset()}
}

\value{ExpressionSet containing the significant GISTIC output    }

\author{Markus Riester}




\examples{
library(copynumbR)
clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
eset <-
copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
    "all_lesions.conf_95.txt", package="copynumbR"), clinical)

# list all recurrent gains and losses
featureData(eset)$Descriptor

# show the copy number distributions of all recurrent alterations
boxplot(t(exprs(eset)))
}
