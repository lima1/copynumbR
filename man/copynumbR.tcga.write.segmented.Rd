\name{copynumbR.tcga.write.segmented}
\alias{copynumbR.tcga.write.segmented}
\title{copynumbR tcga write segmented}
\description{Create a segmented copy number file from TCGA Level 3 data}
\usage{copynumbR.tcga.write.segmented(path = ".", file = "tcga.seg", 
    filter = "\\\.hg18", level = 3, sdrf.file = NULL, sdrf.barcode.col = 1, 
    sdrf.filename.col = "Derived.Array.Data.Matrix.File", verbose = TRUE, 
    ...)}
\arguments{
  \item{path}{The path containing the Level_3 folder.}
  \item{file}{The segemented data output file.}
  \item{filter}{Consider only files passing the specified grep filter. Default the genome version.}
  \item{level}{The TCGA access level.}
  \item{sdrf.file}{Translate filename to TCGA barcode with specified sdrf file.}
  \item{sdrf.barcode.col}{Column with the barcode.}
  \item{sdrf.filename.col}{Column with the data matrix.}
  \item{verbose}{Print some additional progress information.}
  \item{\dots}{Additional arguments passed to read.delim() for reading the files.}
}

\value{A data.frame containing the segmented data.   }

\author{Markus Riester}





