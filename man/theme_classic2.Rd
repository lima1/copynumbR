\name{theme_classic2}
\alias{theme_classic2}
\title{theme classic2}
\description{A simple ggplot2 theme, no fancy stuff.}
\usage{theme_classic2(base_size = 12, base_family = "")}
\arguments{
  \item{base_size}{Theme base font size.}
  \item{base_family}{Theme base font family.}
}



\author{Markus Riester}




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

plot(p+theme_classic2(14)) 
}
