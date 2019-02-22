set.seed(999)
n <- 1000
df <- data.frame(factor = sample(letters[1:8], n, replace = TRUE), 
                 x = rnorm(n), y = runif(n))
circos.par("track.height = 0.1")
circos.initialize(factors = df$factor, x = df$x)
# layer 1
circos.track(factors = df$factor, y = df$y,
             panel.fun = function(x,y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + uy(5, "mm"),
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
col <- rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(df$factor, df$x, df$y, col = col, pch =16,
                   cex = 0.5)
circos.text(-1, 0.5, "haha", sector.index = "a", track.index = 2)
# layer 2
bgcol <- rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(df$factor, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)
# layer 3
circos.track(factors = df$factor, x = df$x, y = df$y,
             panel.fun = function(x, y) {
               ind = sample(length(x), 10)
               x2 = x[ind]
               y2 = y[ind]
               od = order(x2)
               circos.lines(x2[od], y2[od])
             })
# layer update
circos.update(sector.index = "d", track.index = 2,
              bg.col = "#FF8080", bg.border = "black")
circos.points(x = -2:2, y = rep(0.5, 5), col = "white")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "updated", col = "white")
# layer 4
circos.track(ylim = c(0,1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 0.1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = rand_color(n_breaks), border = NA)
})
# layer 5
circos.link("a",0,"b",0,h = 0.4)
circos.link("c", c(-0.5,0.5), "d", c(-0.5,0.5), col = "red",
            border = "blue", h = 0.2)
circos.link("e", 0, "g", c(-1,1), col = "green", border = "black",
            lwd = 2, lty = 2)
circos.link("a",c(-2,-1),"e",c(-2,-1))

circos.initialize(test$chrom,test$position)
circos.track(test$chrom, x = test$position, y = test$hits, ylim = c(1,30), 
             panel.fun = function(x,y) {
               circos.lines(x,y,type = "h",baseline = 20, col = ifelse(y >=20, "red", "blue"))
               circos.text(CELL_META$xcenter,CELL_META$ycenter, CELL_META$sector.index,niceFacing = TRUE)
             })

# genomic track
chr_bins <- tileGenome(genome[c("chr1","chr6")], tilewidth = 20000, cut.last.tile.in.chrom = TRUE)
rep1_chr <- BinChipseq(rep1, chr_bins)
test1 <- data.frame(chr=rep1_chr1@seqnames,start=start(rep1_chr1),end=end(rep1_chr1),value=rep1_chr1$score)
test6 <- data.frame(chr=rep1_chr6@seqnames,start=start(rep1_chr6),end=end(rep1_chr6),value=rep1_chr6$score)
test <- rbind(test1, test6)
circos.initializeWithIdeogram(species = "mm10", chromosome.index = c("chr1"))
circos.genomicTrack(test, numeric.column= 4, bg.border = NA,
                    panel.fun = function(region, value,...) {
                      circos.genomicLines(region,value,..., type = "h")
                    })
circos.genomicLink(test_start, test_end)
circos.genomicTrack(test_np6, numeric.column = 4, bg.border = NA,
					panel.fun = function(region, value, ...) {
						circos.genomicRect(region, value, ..., ytop = CELL_META$ylim[2]-1, ybottom = CELL_META$ylim[1]+1)
					})
