suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
library(Matrix)

## load dmat_all, dataset df (ddf), color vector (col)
source("code/compare/data_helpers.R")

## convert distance matrix into "edge" list
unpack <- function(x) {
  # x => distance matrix
  ex <- Matrix::Matrix(as.matrix(x), sparse=TRUE)
  es <- summary(ex)
  es[,1] <- dtypes[es[,1]]
  es[,2] <- dtypes[es[,2]]
  ## filter identical data types
  es <- es[es[,1]==es[,2],]
  colnames(es) <- c("dtype", "dtype", "dist")
  as.data.frame(es[,2:3])
}

n <- c(25, 50, 100, 200, 350, 700, 1800, 4300, 9000)
dist_plot <- function(X, metadf, main="") {
  Xdf <- left_join(unpack(X), metadf, by='dtype')
  ggplot(aes(x=n, y=dist, col=method), data=Xdf) +
    stat_summary(fun.y=mean, geom='line') +
    stat_summary(fun.data=mean_sdl, geom='linerange') +
    scale_color_manual(values=col) +
    scale_x_log10(breaks=n) +
    ylab("distance") + ggtitle(main) +
    theme_bw()
}


pdf("plots/intramethod_distances.pdf", width=7, height=4)
for (i in 1:length(dmat_all))
 print(dist_plot(dmat_all[[i]], ddf, names(dmat_all)[i]))
dev.off()
