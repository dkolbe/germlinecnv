#!/usr/bin/env Rscript

# Environment setup
library(AnnotationDbi)

# Template-interpolated by main.nf
buildID <- "${ucsc_genome}"

# Error checking on parameters

# download and save based on buildID

txdb = paste("txdb.refseq.",buildID,".sqlite")

ucsc_refseq <- makeTxDbFromUCSC(genome=buildID,tablename="refGene")
saveDb(ucsc_refseq,txdb)

