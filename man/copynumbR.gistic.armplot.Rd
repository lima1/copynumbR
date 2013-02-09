\name{copynumbR.gistic.armplot}
\alias{copynumbR.gistic.armplot}
\title{copynumbR gistic armplot}
\description{Plot chromosome numbers as estimated by GISTIC}
\usage{copynumbR.gistic.armplot(file = "broad_values_by_arm.txt", subtypes = NULL, 
    chromosomes = c(1:12, 16:20), ncol = NULL, ...)}
\arguments{
  \item{file}{The GISTIC output file broad_values_by_arm.txt}
  \item{subtypes}{If not NULL, then stratify by subtype}
  \item{chromosomes}{Display only the specified chromosomes. }
  \item{ncol}{Number of columns in the plot}
  \item{\dots}{
}
}

\value{A ggplot2 object}

\author{Markus Riester}




\examples{
library(copynumbR)

p <-
    copynumbR.gistic.armplot(file=
        system.file("extdata/gistic_stransky_bladder", 
            "broad_values_by_arm.txt", package="copynumbR")
    )

plot(p)
}
