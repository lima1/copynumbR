#!/bin/sh

R --vanilla <<RSCRIPT
library(inlinedocs);
package.skeleton.dx(".")
RSCRIPT

cd ..
R CMD install --build copynumbR
cd copynumbR
