# Environment setup
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)

# FIXME strongly consider making single run per sample rather than whole batch

# FIXME set to the work directory
setwd("~/Desktop/temp/20220112-NGS687/")

# Load standard tracks and annotation. We are retrieving UCSC refseq, so can be slow and may update on subsequent runs.
gtrack <- GenomeAxisTrack(showTitle=TRUE,rotation.title=0,cex=1.2)

# FIXME parameterize genome options
if(file.exists('txdb.refseqhg19.sqlite')) {
  ucsc_refseq = loadDb(file = 'txdb.refseqhg19.sqlite')
} else {
  ucsc_refseq <- makeTxDbFromUCSC(genome="hg19",tablename="refGene")
  saveDb(ucsc_refseq,'txdb.refseqhg19.sqlite')
}  

# FIXME parameterize region list input
probelist <- read.table("otoscope_v9_probes_CNVsubset.bed")
# FIXME make fewer assumptions about columns
names(probelist) <- c("CHR","START","END","GENE","INHERITANCE","LOCUS","TRANSCRIPT","TIER")
probelist_gr <- makeGRangesFromDataFrame(probelist,seqnames.field="CHR",start.field="START",end.field="END",keep.extra.columns = TRUE)
# Add gene marker regions to problist
# FIXME do this in a configurable, not hard-coded way
probelist_gr <- append(probelist_gr,GRanges("15",IRanges(start=43891678,end=44010554),GENE="STRCregion"))
probelist_gr <- append(probelist_gr,GRanges("16",IRanges(start=21675224,end=22588189),GENE="OTOAregion"))
probelist_gr <- append(probelist_gr,GRanges("13",IRanges(start=20761532,end=21106120),GENE="GJB2region"))
probelist_gr <- append(probelist_gr,GRanges("22",IRanges(start=19136435,end=21356679),GENE="TBX1region"))

# FIXME
genome(seqinfo(probelist_gr)) <- "GRCh37"

# Load aggregated copy-number data for a batch
# FIXME parameterize input
agg <- read.table("aggregated_copy_numbers.txt",sep="\t",header=TRUE)
agg_gr <- makeGRangesFromDataFrame(agg,seqnames.field="CONTIG",start.field="START",end.field="END",keep.extra.columns = TRUE)
# FIXME 
genome(seqinfo(agg_gr)) <- "GRCh37"

# FIXME
ids <- read.table("aggregated_copy_numbers.txt",sep="\t",header=FALSE,nrows=1,as.is=TRUE)[-(1:3)]

# Load segmented call data
for (samplename in ids) {
  # FIXME match updated file naming pattern
  segments <- read.table(paste("genotyped-segments-",samplename,".vcf",sep=""),comment.char = "#")
  # FIXME
  names(segments) <- c("CHR","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","DATA")
  segments <- segments[which(segments$ALT != '.'),]
  if (isEmpty(segments)) { next; }
  segments$END <- as.numeric(gsub("END=","",segments$INFO))
  segments_gr <- makeGRangesFromDataFrame(segments,seqnames.field="CHR",
                                          start.field="POS",end.field="END",
                                          keep.extra.columns = TRUE)
  # FIXME
  genome(seqinfo(segments_gr)) <- "GRCh37"
  
  genelist <- unique(mcols(subsetByOverlaps(probelist_gr,segments_gr))$GENE)
  
  # Expand gene list to cover regions
  # FIXME do this in a configurable, not hard-coded way
  if (any(c('STRC','STRCP1') %in% genelist)) {
    genelist <- union(genelist,c('STRC','CATSPER2','STRCP1','STRCregion'))
  }
  if (any(c('OTOA') %in% genelist)) {
    genelist <- union(genelist,c('OTOA','OTOAP1','OTOAregion'))
  }
  if (any(c('GJB2','GJB6','CRYL1') %in% genelist)) {
    genelist <- union(genelist,c('GJB2','GJB6','CRYL1','GJB2region'))
  }
  if (any(c('TBX1') %in% genelist)) {
    genelist <- union(genelist,c('TBX1','GSC2','RTN4R','ZNF74','THAP7','TBX1region'))
  }
  
  for (gene in genelist) {
    if (is.na(gene)) { next; } # Might need to deal with these, but later
    extent <- probelist_gr[which(probelist_gr$GENE==gene),]+50
  
    backgroundtrack <- DataTrack(subsetByOverlaps(agg_gr,extent),
                                 name="Copies",col="gray",ylim=c(0,5),
                                 baseline=c(0,1,2,3,4,5),cex=.5)
    current <- agg_gr[,make.names(samplename)]
    current <- subsetByOverlaps(current,extent)
    sampletrack <- DataTrack(current,col="black",ylim=c(0,5),cex=1)
    deltrack <- DataTrack(subsetByOverlaps(current,
                                           segments_gr[which(segments_gr$ALT == "<DEL>"),]),
                        col="red",ylim=c(0,5),cex=1,legend=FALSE)
    duptrack <- DataTrack(subsetByOverlaps(current,
                                           segments_gr[which(segments_gr$ALT == "<DUP>"),]),
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
    # FIXME check output name
    png(file=paste(samplename,gene,".png",sep=""),width=1200,height=900)
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
