circos.clear()
circos.genomicInitialize(rna2)
circos.genomicTrackPlotRegion(rna2, ylim = c(0,1), track.height = 0.05, bg.border = 'white',
                              panel.fun = function(region, value, ...) {
                                for (i in (1: nrow(region))) {
                                  reg <- region[i,]
                                  val <- value[1,]
                                  color <- ifelse(i %% 2 == 1, 'black', 'grey')
                                  circos.genomicRect(reg, val, col = color, border = 'white')
                                }
                              })
circos.genomicLink(links_from, links_to, col = rand_color(nrow(links_from), transparency = 0.7))

# arrange
re_genome <- function(df, bin){
  end <- vector()
  for (i in range(1, nrow(df))) {
    end <- append(end, (df$end[i] - df$start[i]) %/% bin + 1)
  }
  df$start <- rep(1, nrow(df))
  df$end <- end
  return(df)
}


# arrange start & end as bin position
re_position <- function(df, start, bin){
  start_bin <- vector()
  end_bin <- vector()
  for (i in range(1, nrow(df))) {
    start_bin <- append(start_bin, (df$start[i] - start) %/% bin + 1)
    end_bin <- append(end_bin, (df$end[i] - start) %/% bin + 1)
  }
  df$start <- start_bin
  df$end <- end_bin
  return(df)
}
