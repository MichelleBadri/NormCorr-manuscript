suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
library(Matrix)

## load dmat_all, dataset df (ddf), color vector (col)
source("code/compare/data_helpers.R")

## convert distance matrix into "edge" list
unpack_cor <- function(x) {
  # x => distance matrix
  ex <- Matrix::Matrix(as.matrix(x), sparse=TRUE)
  es <- summary(ex)
  es[,1] <- dtypes_cor[es[,1]]
  es[,2] <- dtypes_cor[es[,2]]
  ## filter identical data types
  es <- es[es[,1]==es[,2],]
  colnames(es) <- c("dtype", "dtype", "dist")
  as.data.frame(es[,2:3])
}

n <- c(25, 50, 100, 200, 350, 700, 1800, 4300, 9000)
dist_plot_cor <- function(X, metadf, main="") {
  Xdf <- left_join(unpack_cor(X), metadf, by='dtype')
   ggplot(aes(x=n, y=dist, col=method, group=group, linetype=association, shape=method), data=Xdf) +
    stat_summary(fun.y=mean,  geom='point', data= subset(Xdf,association!="cor"), size=2.5) +
    stat_summary(fun.y=mean,  geom='line') +
    stat_summary(fun.data=mean_sdl, geom='linerange', alpha=0.4) +
   # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
    scale_color_manual(values=col) +
    scale_shape_manual(values=as.numeric(shape)) +
    scale_linetype_manual(values=c("dashed","solid","twodash"))+
    scale_x_log10(breaks=n) +
    ylab("distance") + ggtitle(main) +
    theme_bw()
}

pdf("plots/intramethod_distances_pearson_corshrink_propr.pdf", width=7, height=5)
for (i in 1: length(dmat_all_cor))
 print(dist_plot_cor(dmat_all_cor[[i]], ddf_cor, names(dmat_all_cor)[i]))
dev.off()



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
    ggplot(aes(x=n, y=dist, col=method, shape=method), data=Xdf) +
      stat_summary(fun.y=mean,  geom='point', size=2.5) +
      stat_summary(fun.y=mean,  geom='line') +
      stat_summary(fun.data=mean_sdl, geom='linerange', alpha=0.4) +
      # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
      scale_color_manual(values=col) +
      scale_shape_manual(values=as.numeric(shape)) +
      scale_x_log10(breaks=n) +
      ylab("distance") + ggtitle(main) +
      theme_bw()
  }

pdf("plots/intramethod_distances.pdf", width=7, height=4)
for (i in 1:length(dmat_all))
  print(dist_plot(dmat_all[[i]], ddf, names(dmat_all)[i]))
dev.off()
