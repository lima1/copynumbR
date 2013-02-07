\name{copynumbR.gistic.focal}
\alias{copynumbR.gistic.focal}
\title{copynumbR gistic focal}
\description{Create a data.frame with all GISTIC focal alterations. Useful for
presentation in knitR/Sweave.}
\usage{copynumbR.gistic.focal(eset, gistic.lesions.file.amp = "amp_genes.conf_95.txt", 
    gistic.lesions.file.del = "del_genes.conf_95.txt", gain = 0.1, 
    loss = -0.1)}
\arguments{
  \item{eset}{The GISTIC lesions file read with copynumbR.gistic.read.lesions}
  \item{gistic.lesions.file.amp}{The GISTIC output file amp_genes.conf_95.txt}
  \item{gistic.lesions.file.del}{The GISTIC output file del_genes.conf_95.txt}
  \item{gain}{Minimum log2 ratio for a copy number gain}
  \item{loss}{Maximum log2 ratio for a copy number loss}
}



\author{Markus Riester}




\examples{
library(copynumbR)
clinical <- read.csv(system.file("extdata", "stransky_bladder_clinical.csv", package="copynumbR"))
eset <-
copynumbR.gistic.read.lesions(system.file("extdata/gistic_stransky_bladder",
    "all_lesions.conf_95.txt", package="copynumbR"), clinical)

focal <-
    copynumbR.gistic.focal(eset, 
        system.file("extdata/gistic_stransky_bladder", 
            "amp_genes.conf_95.txt", package="copynumbR"),
        system.file("extdata/gistic_stransky_bladder",
            "del_genes.conf_95.txt", package="copynumbR")
    )
}
