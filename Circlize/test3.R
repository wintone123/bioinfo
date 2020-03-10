arrange_genes <- function(df) {
    df_from <- df %>% group_by(gene1, position1) %>% summarise(counts = sum(counts))
    df_to <- df %>% group_by(gene2, position2) %>% summarise(counts = sum(counts))
    df_out <- data.frame(gene = append(df_from$gene1, df_to$gene2),
                         from = as.numeric(append(df_from$position1, df_to$position2)),
                         counts = as.numeric(append(df_from$counts, df_to$counts)),
                         stringsAsFactors = FALSE)
    df_out$to <- df_out$from + 100
    df_out <- select(df_out, c(1,2,4,3))
    return(df_out)
}


arrange_links <- function(df) {
    df_out <- select(df, c(1,2,5)) 
    df_out$from1 <- as.numeric(df$position1) + sample(100, nrow(df), replace = TRUE)
    df_out$to1 <- df_out$from1
    df_out$from2 <- as.numeric(df$position2) + sample(100, nrow(df), replace = TRUE)
    df_out$to2 <- df_out$from2
    df_out <- select(df_out, c(1, 4, 5, 2, 6, 7, 3))
    return(df_out)
}


interaction_graph <- function(genome, links) {
  circos.clear()
  circos.genomicInitialize(genome)
  # counts show as histograms
  circos.genomicTrackPlotRegion(genome, ylim = c(0,1), track.height = 0.1, bg.border = '#ffffffff',
                                panel.fun = function(region, value, ...) {
                                  for (i in seq_len(nrow(region))) {
                                    reg <- region[i,]
                                    val <- value[i,] / max(value[1])
                                    circos.genomicRect(reg, val, ybottom = 0, ytop = val, border = '#00000090')
                                  }
                                })
  # genome tracks
  circos.genomicTrackPlotRegion(genome, ylim = c(0,1), track.height = 0.05, bg.border = 'white',
                                panel.fun = function(region, value, ...) {
                                  for (i in seq_len(nrow(region))) {
                                    reg <- region[i,]
                                    val <- value[i,]
                                    color <- ifelse(i %% 2 == 1, '#000000ff', '#00000080')
                                    circos.genomicRect(reg, val, col = color, border = '#ffffffff')
                                  }
                                })
  # genome links
  circos.genomicLink(links[1:3], links[4:6], 
                     col = ifelse(links$counts / sum(links$counts) >= 0.01, '#ff0077de', '#008cffce'), 
                     lwd = 0.7)
}
