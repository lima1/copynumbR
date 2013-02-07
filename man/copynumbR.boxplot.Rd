\name{copynumbR.boxplot}
\alias{copynumbR.boxplot}
\title{copynumbR boxplot}
\description{A boxplot showing correlation of copy number and Expression for matched
data}
\usage{copynumbR.boxplot(eset.cn, eset.expr, cutoffs = c(-Inf, -1.3, 
    -0.1, 0.1, 0.9, Inf), cutoff.labels = c("High Loss", "Loss", 
    "Normal", "Gain", "Amplification"), probesets = NULL, min.samples = 3, 
    sqrt = FALSE, highlight = NULL, highlight.labels = NULL, 
    xlab = "Copy Number", ylab = "Expression")}
\arguments{
  \item{eset.cn}{ExpressionSet with copy number data}
  \item{eset.expr}{ExpressionSet with expression data}
  \item{cutoffs}{Copy number cutoffs }
  \item{cutoff.labels}{The labels of these cutoffs}
  \item{probesets}{Display only these genes. If null, show all genes.}
  \item{min.samples}{Minimum number of samples in each cutoff category}
  \item{sqrt}{Square root transform the data for read counts?}
  \item{highlight}{highlight some samples (not yet implemented)}
  \item{highlight.labels}{The label of these highlighted samples shown in the legend}
  \item{xlab}{The label of the x-axis}
  \item{ylab}{The label of the y-axis}
}

\value{A ggplot2 object.    }

\author{Markus Riester}





