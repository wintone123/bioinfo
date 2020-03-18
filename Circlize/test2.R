interaction_graph <- function(genome, from, to) {
  circos.clear()
  circos.genomicInitialize(genome)
  # counts show as histograms
  circos.genomicTrackPlotRegion(genome, ylim = c(0,1), track.height = 0.1, bg.border = 'white',
                                panel.fun = function(region, value, ...) {
                                  for (i in seq_len(nrow(region))) {
                                    reg <- region[i,]
                                    val <- value[i,] / max(value[1])
                                    circos.genomicRect(reg, val, ybottom = 0, ytop = val, border = 'white', 
                                                       col = ifelse(val >= 0.5, 'red', 'green'))
                                  }
                                })
  # genome tracks
  circos.genomicTrackPlotRegion(genome, ylim = c(0,1), track.height = 0.05, bg.border = 'white',
                                panel.fun = function(region, value, ...) {
                                  for (i in seq_len(nrow(region))) {
                                    reg <- region[i,]
                                    val <- value[i,]
                                    color <- ifelse(i %% 2 == 1, 'black', 'grey')
                                    circos.genomicRect(reg, val, col = color, border = 'white')
                                  }
                                })
  # genome links
  circos.genomicLink(from, to)
}


interaction_graph <- function(genome, links) {
  circos.clear()
  circos.genomicInitialize(genome)
  # counts show as vertical lines
  circos.genomicTrackPlotRegion(genome, ylim = c(0,1), track.height = 0.1, bg.border = 'white',
                                panel.fun = function(region, value, ...) {
                                  for (i in seq_len(nrow(region))) {
                                    reg <- region[i,]
                                    val <- value[i,] / max(value[1])
                                    circos.genomicLines(reg, val, numeric.column = 1, type = 'h', lwd = 1,
                                                        col = ifelse(val >= 0.5, 'red', 'green'))
                                  }
                                })
  # genome tracks
  circos.genomicTrackPlotRegion(genome, ylim = c(0,1), track.height = 0.05, bg.border = 'white',
                                panel.fun = function(region, value, ...) {
                                  for (i in seq_len(nrow(region))) {
                                    reg <- region[i,]
                                    val <- value[i,]
                                    color <- ifelse(i %% 2 == 1, '#FF00FF60', 'grey')
                                    circos.genomicRect(reg, val, col = color, border = NA)
                                  }
                                })
  # genome links
  circos.genomicLink(links[1:3], links[4:6])
}


# arrange
re_genome <- function(df, bin) {
  end <- vector()
  for (i in seq_len(nrow(df))) {
    end <- append(end, (df$end[i] - df$start[i]) %/% bin + 1)
  }
  df$start <- rep(1, nrow(df))
  df$end <- end
  return(df)
}


# arrange start & end as bin position
re_position <- function(df, start, bin) {
  start_bin <- vector()
  end_bin <- vector()
  for (i in seq_len(nrow(df))) {
    start_bin <- append(start_bin, (df$start[i] - start) %/% bin + 1)
    end_bin <- append(end_bin, (df$end[i] - start) %/% bin + 1)
  }
  df$start <- start_bin
  df$end <- end_bin
  return(df)
}
