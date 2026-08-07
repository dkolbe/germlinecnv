#!/usr/bin/env Rscript

# Environment setup
suppressPackageStartupMessages({
library(GenomicFeatures)
library(AnnotationDbi)
library(txdbmaker)
})

# Template-interpolated by main.nf
buildID <- "${ucsc_genome}"
txdb_path <- "${txdb_name}"

# Error checking on parameters

# download and save based on buildID

ucsc_refseq <- txdbmaker::makeTxDbFromUCSC(genome=buildID,tablename="refGene")
saveDb(ucsc_refseq,txdb_path)

