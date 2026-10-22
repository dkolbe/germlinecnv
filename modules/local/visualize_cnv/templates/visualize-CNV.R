#!/usr/bin/env Rscript

# Environment setup
suppressPackageStartupMessages({
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)
})

# FIXME strongly consider making single run per sample rather than whole batch

# Template-interpolated by main.nf
calls      <- "${calls_table}"
counts     <- "${counts}"
bg_counts  <- "${aggregated_counts}"
regions    <- "${regions}"
txdb       <- "${txdb}"
buildID    <- "${ucsc_genome}"
graph_all  <- "${graph_all}"
graph_some  <- "${graph_some}"

graph_some_list <- unlist(strsplit(graph_some,','))

# Error checking on parameters
sprintf("Calls: %s",calls)
sprintf("Counts: %s",counts)
sprintf("Background counts: %s",bg_counts)
sprintf("Regions: %s",regions)
sprintf("txdb: %s",txdb)
sprintf("Build: %s",buildID)
sprintf("graph_all: %s",graph_all)
sprintf("graph_some: %s",graph_some)

# Load standard tracks and annotation. We are retrieving UCSC refseq, so can be slow and may update on subsequent runs.
gtrack <- GenomeAxisTrack(showTitle=TRUE,rotation.title=0,cex=1.2)

if(file.exists(txdb)) {
  ucsc_refseq = loadDb(file = txdb)
} else {
  # ERROR
  warning("nf-core module implementation requires txdb to be set up in advance, use DOWNLOAD_DB")
  quit(status=1)
}  

probelist <- read.table(regions)
probelist <- probelist[,1:4]
names(probelist) <- c("CHR","START","END","GENE")
probelist_gr <- makeGRangesFromDataFrame(probelist,seqnames.field="CHR",start.field="START",end.field="END",keep.extra.columns = TRUE)
genome(seqinfo(probelist_gr)) <- buildID

# Load aggregated copy-number data for a batch
agg <- read.table(bg_counts,sep="\t",header=TRUE)
agg_gr <- makeGRangesFromDataFrame(agg,seqnames.field="CONTIG",start.field="START",end.field="END",keep.extra.columns = TRUE)
genome(seqinfo(agg_gr)) <- buildID

# Also build the sample-specific counts into a GR object
sample_counts <- read.table(counts,sep="\t",comment.char = "@",header=TRUE)
sample_counts_gr <- makeGRangesFromDataFrame(sample_counts,seqnames.field="CONTIG",start.field="START",end.field="END",keep.extra.columns = TRUE)
genome(seqinfo(sample_counts_gr)) <- buildID

# FIXME should be able to just pull this from the table read above rather than reading the file again
# ... but reading just the first line is NBD
# ids <- read.table(bg_counts,sep="\t",header=FALSE,nrows=1,as.is=TRUE)[-(1:3)]
# Replaced with only processing the single sample passed with metadata rather than all in the aggregated counts dataset

# Load segmented call data
# for (samplename in ids) {
for (samplename in "${meta.id}") {
  segments <- read.table(calls,sep="\t",header=TRUE)
  if (isEmpty(segments)) { next; }
  segments_gr <- makeGRangesFromDataFrame(segments,seqnames.field="CONTIG",
                                          start.field="START",end.field="END",
                                          keep.extra.columns = TRUE)
  genome(seqinfo(segments_gr)) <- buildID
  
  if (graph_all) {
    genelist2 <- unique(mcols(probelist_gr)\$GENE)
  } else {
    genelist1 <- unique(mcols(subsetByOverlaps(probelist_gr,segments_gr))\$GENE)

    # Exapnd gene list to include specific requested genes
    genelist1 <- union(genelist1,graph_some_list)

    # Expand gene list to cover all genes within overview regions
    genelist_regions_gr <- probelist_gr[elementMetadata(probelist_gr)[,"GENE"] %in% genelist1]
    genelist2 <- unique(mcols(subsetByOverlaps(probelist_gr,genelist_regions_gr))\$GENE)
  }
  
  for (gene in genelist2) {
    if (is.na(gene)) { next; } # Might need to deal with these, but later
    extent <- probelist_gr[which(probelist_gr\$GENE==gene),]+50
  
    backgroundtrack <- DataTrack(subsetByOverlaps(agg_gr,extent),
                                 name="Copies",col="gray",ylim=c(0,5),
                                 baseline=2,cex=.5)
    #current <- agg_gr[,make.names(samplename)]
    current <- subsetByOverlaps(sample_counts_gr,extent)
    sampletrack <- DataTrack(current,col="black",ylim=c(0,5),cex=1)
    deltrack <- DataTrack(subsetByOverlaps(current,
                                           segments_gr[which(segments_gr\$ALT == "<DEL>"),]),
                        col="red",ylim=c(0,5),cex=1,legend=FALSE)
    duptrack <- DataTrack(subsetByOverlaps(current,
                                           segments_gr[which(segments_gr\$ALT == "<DUP>"),]),
                        col="blue",ylim=c(0,5),cex=1,legend=FALSE)
    overtrack <- OverlayTrack(list(backgroundtrack,sampletrack,deltrack,duptrack))
    names(gtrack) <- as.vector(runValue(seqnames(extent)))
    chromosome(gtrack) <- as.vector(runValue(seqnames(extent)))
    txTr <- GeneRegionTrack(ucsc_refseq, name="Genes",
                            chromosome = as.vector(runValue(seqnames(extent))),
                            start=min(start(extent)),
                            end=max(end(extent)),
                            cex=2)
    if (length(unique(transcript(txTr)))>10) { stacking(txTr) <- "dense" }
    png(file=paste(samplename,'.',gene,".png",sep=""),width=1200,height=900)
    plotTracks(list(gtrack,overtrack,txTr),
               main=paste(samplename,gene),
               from=min(start(extent)),to=max(end(extent)),
               chromosome = as.vector(runValue(seqnames(extent))),
               transcriptAnnotation="symbol",just.group="below",
               col.title="black",col.axis="black",
               title.width=1.5,cex.title=1.5,cex.axis=1,sizes=c(1,5,2))
    dev.off()
  }
}
