set.seed(2333)
mat <- matrix(sample(18, 18), 3, 6)
rownames(mat) <- paste0('E', 1 : 3)
colnames(mat) <- paste0('S', 1 : 6)

df <- data.frame(from = rep(rownames(mat), times = ncol(mat)), 
                 to = rep(colnames(mat), each = nrow(mat)),
                 value = as.vector(mat),
                 stringsAsFactors = FALSE)


# output matrix
df_arrange <- function(df) {
    mat <- matrix(0, length(unique(df$type1)), length(unique(df$type2)))
    rownames(mat) <- unique(df$type1)
    colnames(mat) <- unique(df$type2)
    for (i in seq_len(nrow(df))) {
        mat[df[i, 1], df[i, 3]] <- mat[df[i, 1], df[i, 3]] + 1
    }
    return(mat)
}


# output dataframe
df_arrange <- function(df) {
    mat <- matrix(0, length(unique(df$type1)), length(unique(df$type2)))
    rownames(mat) <- unique(df$type1)
    colnames(mat) <- unique(df$type2)
    for (i in seq_len(nrow(df))) {
        mat[df[i, 1], df[i, 3]] <- mat[df[i, 1], df[i, 3]] + 1
    }
    df <- data.frame(from = rep(rownames(mat), times = ncol(mat)),
                    to = rep(colnames(mat), each = nrow(mat)),
                    value = as.vector(mat),
                    stringsAsFactors = FALSE)
    return(df)
}


# gather reverse items
re_remove <- function(df) {
    df <- arrange(df, from)
    a <- 0
    while (TRUE) {
        if (a == 0) {
            df_out <- df[1,]
            a <- 1
        } else {
            df_out <- rbind(df_out, df[1,])
        }
        df <- df[-1,]
        if (nrow(df) == 0) {
            break
        }
        n <- nrow(df_out)
        for (i in seq_len(nrow(df))) {
            if (df$to[i] == df_out$from[n] & df$from[i] == df_out$to[n]) {
                df_out$value[n] <- df_out$value[n] + df$value[i]
                df <- df[-i,]
                break
            }
        }
        if (nrow(df) == 0) {
            break
        }
    }
    return(df_out)
}

chordDiagram()