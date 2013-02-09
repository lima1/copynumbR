\name{copynumbR.eset}
\alias{copynumbR.eset}
\title{copynumbR eset}
\description{Generate an ExpressionSet from a data.frame containing the data and a
data.frame containing clinical data.}
\usage{copynumbR.eset(data, clinical, data.col = 2, id.col = 1, post.process.fun = NULL, 
    ...)}
\arguments{
  \item{data}{A data.frame containing genomic data (expression, copy number, ...).}
  \item{clinical}{A data.frame containing matched clinical data.}
  \item{data.col}{Sample data starts at this column.}
  \item{id.col}{The column of the matching sample ids in the clinical data.frame.}
  \item{post.process.fun}{A function to manipulate the final ExpressionSet (adding a Surv obj for
example).}
  \item{\dots}{Additional arguments passed to post.process.fun().}
}



\author{Markus Riester}




\examples{
library(copynumbR)
clinical <- read.csv(system.file("extdata", 
    "stransky_bladder_clinical.csv", package="copynumbR"))

data <- read.delim(system.file("extdata/gistic_stransky_bladder",
"broad_values_by_arm.txt", package="copynumbR"),stringsAsFactors=FALSE)

# add an example post.process.fun. makes it easy to change the annotation
.curateGender <- function(eset, ...) {
    eset$GENDER.2 <- as.factor(as.character(ifelse(eset$GENDER=="M", "male",
        "female")))
    eset
}

eset <- copynumbR.eset(data, clinical, data.col=2, post.process.fun=.curateGender)

eset$GENDER.2
}
