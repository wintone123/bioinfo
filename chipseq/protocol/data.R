library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(biomaRt)

# read extension
prepareCHIPseq <- function(reads){
  frag.len = median( estimate.mean.fraglen(reads) )
  cat( paste0( 'median size is ', round(frag.len)))
  reads.extended = resize(reads, width = frag.len)
  return( trim(reads.extended) )
}

# bins generation
binsize <- 200
bins <- tileGenome(c(chr6 = 149517037), tilewidth = binsize, cut.last.tile.in.chrom = TRUE)

# bining
BinChipseq <- function(reads,bins){
  mcols(bins)$score=countOverlaps(bins,reads)
  return(bins)
}

input_bins <- BinChipseq(input, bins)
rep1_bins <- BinChipseq(rep1, bins)
rep2_bins <- BinChipseq(rep2, bins)

#binned data export for IGV
export(input_bins, con = "input_chr6_bedgraph", format = "bedgraph")

# obtaining objects si from mm9
library(BSgenome.Mmusculus.UCSC.mm9)
genome <- BSgenome.Mmusculus.UCSC.mm9
si <- seqinfo(genome)
si <- si[paste0("chr", c(1:19, "x", "y"))]

# obtaining objects bm from mm9
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",
                host = "jul2018.archive.ensembl.org")
fm <- Gviz:::.getBMFeatureMap()
fm["symbol"] <- "external_gene_id"
bm <- BiomartGeneRegionTrack(chromosome = "chr6", genome = "mm9",
                             start = 34000000, end = 35000000,
                             biomart = mart, filter = list("with_ox_refseq_mrna"=TRUE),
                             size = 4, name = "Refseq", utr5 = "red3", utr3 = "red3",
                             protein_coding = "black", col.line = NULL, cex = 7,
                             collapseTranscript = "longest", featureMap = fm)
AT <- GenomeAxisTrack()

# plot bm 
plotTracks(c(bm,AT), from = 34000000, to = 35000000,
           transcriptAnnotation = "symbol", window = "auto",
           cex.title = 1, fontsize = 10)

input_track <- DataTrack(input_bins, strand = "*", genome = "mm9", 
                                        col.histogram = "gray", fill.histrogram = "black",
                                        name = "input", col.axis = "black", cex.axis = 0.4, ylim = c(0,150))
rep1_track <- DataTrack(rep1_bins, strand = "*", genome = "mm9", 
                                        col.histogram = "gray", fill.histrogram = "black",
                                        name = "input", col.axis = "black", cex.axis = 0.4, ylim = c(0,150))
rep2_track <- DataTrack(rep2_bins, strand = "*", genome = "mm9", 
                                        col.histogram = "gray", fill.histrogram = "black",
                                        name = "input", col.axis = "black", cex.axis = 0.4, ylim = c(0,150))

# plot tracks with genomic features
plotTracks(c(input_track, rep1_track, rep2_track, bm,AT), from = 122630000, 
           to = 122700000, transcriptAnnotation = "symbol",window = "auto", 
           type = "histogram", cex.title = 0.7,fontsize = 10)

# CHIP seq peaks
test_peaks <- import.bed(file.path("test_summits.bed"))
test_peaks_track <- AnnotationTrack(test_peaks, genome = "mm9", name = "test peak", chromosome = "chr6",
                              shape = "box", fill = "blue3", size = 2)
plotTracks(c(test_peaks_track, bm, AT), from = 34000000, to = 35000000, 
           transcriptAnnotation = "symbol", window = "auto",
           type = "histogram", cex.title = 0.7, fontsize = 10)

# enrichment
enriched_regions <- Reduce(subsetByOverlaps, list(test_peaks))
enr_reg_track <- AnnotationTrack(enriched_regions,
                                 genome = "mm9", name = "Enriched regions",
                                 chromosome = "chr6", shape = "box", fill = "green3", size = 2)

# promoter isolation （mouse)
listAttributes(mart)[1:3,]
ds <- useDataset("mmusculus_gene_ensembl", mart = mart)
chroms <- 6
egs <- getBM(attributes = c("ensembl_gene_id", "external_gene_id",
                            "chromosome_name", "start_position",
                            "end_position", "strand"), 
             filters = "chromosome_name", values = chroms, mart = ds)
egs$TSS <- ifelse(egs$strand == "1", egs$start_position, egs$end_position)

# set +/-200 bp around TSS as promoter region
promoter_regions <- GRange(sequences = Rle(paste0("chr" egs$chromosome_name)),
                           ranges = IRanges(start = egs$TSS - 200,
                                             end = egs$TSS + 200),
                           strand = Rle(rep("*", nrow(egs))),
                           gene = egs$external_gene_id)

# overlapping promoter with enriched regions
ovlp2 <- findOverlaps(enriched_regions, promoter_regions)
cat(sprintf("%d of %d promoters are overlapped by an enriched region.",
            length(unique(ovlp2@subjectHits)), length(promoter_regions)))
ovlp2b <- findOverlaps(promoter_regions, enriched_regions)
cat(sprintf("d% od d% enriched regions overlap a promoter.",
            length(unique(ovlp2b@subjectHits)), length(enriched_regions)))
promoter_total_length <- sum(width(reduac(promoter_regions)))

# which promoter overlapped with peak
pos_TSS <- egs[unique(findOverlaps(promoter_regions, enriched_regions)@queryHits),]

# distribution of peaks around a subset of promoters
tiles <- sapply(1:nrow(pos_TSS), function(i)
    if(pos_TSS$strand[i] == "1")
        pos_TSS$TSS[i] + seq(-1000, 900, length.out = 20)
    else 
        pos_TSS$TSS[i] + seq(900, -1000, length.out = 20))
tiles <- GRanges(tilename = paste(rep(pos_TSS$ensembl_gene_id, each = 20), 1:20, sep = "-"),
                 seqname = Rle(rep(paste0("chr", pos_TSS$chromosome_name), each = 20)),
                 ranges = IRanges(start = as.vector(tiles),
                                  width = 100,)
                 strand = Rle(rep("*", length(as.vector(tiles)))),
                 seqinfo = si)

# count reads mapping to each tile
H3K27ac_p <- countOverlaps(tile, rep1) + countOverlaps(tile, rep2)
H3K27ac_p_matrix <- matrix(H3K27ac_p, nrow = nrow(pos_TSS), ncol = 20, byrow = TRUE)

# output heatmap and plot
colors <- colorRampPalette(c("white", "red", "gray", "black"))(100)
layout(mat = matrix(c(1, 2, 0, 3), 2, 2),
       widths = c(2, 2, 2),
       height = c(0.5, 5, 0.5, 5), TRUE)
par(mar = c(4, 4, 1, 5, 1))
image(seq(0, max(H3K27ac_p_matrix), length.out = 100), 1
      matrix(seq(0, max(H3K27ac_p_matrix), length.out = 100), 100, 1),
      col = colors,
      xlab = "Distance from TSS", ylab = " ",
      main = "Number of reads", yaxt = "n",
      lwd = 3, axes = TRUE)
box(col = "black", lwd = 2)
image(x = seq(-1000, 1000, length.out = 20),
      y = 1:nrow(H3K27ac_p_matrix),
      z = t(H3K27ac_p_matrix[order(rowSums(H3K27ac_p_matrix)),]),
      col = colors,
      xlab = "Distance from TSS (bp)",
      ylab = "Promoters"m lwd = 2)
box(col = "black", lwd = 2)
abline(v = 0, lwd = 1, col = "gray")
plot(x = seq(-1000, 1000, length.out = 20),
     y = conMeans(H3K27ac_p_matrix),
     ty = "b", pch = 19, col = "red4", lwd = 2,
     xlab = "Distance from TSS (bp)",
     ylab = "Mean tag count")
abline(h = seq(1, 100, by = 5),
       v = seq(-1000, 1000, length.out = 20),
       lwd = 0.25, col = "gray")
box(col = "black", lwd = 2)

# biomart release 93
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",
                host = "jul2018.archive.ensembl.org")

listAttributes(mart)[c(1,3,5),]
ds <- useDataset("mmusculus_gene_ensembl", mart = mart)
egs <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                            "chromosome_name", "start_position",
                            "end_position", "strand"), mart = ds)