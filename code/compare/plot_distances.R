suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
library(Matrix)

## load dmat_all, dataset df (ddf), color vector (col)
source("code/compare/data_helpers.R")

n <- yaml::yaml.load_file('code/helpers/subsamples.yml')
dist_plot <- function(X, metadf, main="") {
    Xdf <- left_join(unpack(X, dtypes_all), metadf, by='dtype')
    ggplot(aes(x=n, y=dist,col=factor(method, levels=allord),
               shape=factor(method, levels=allord)), data=Xdf) +
      stat_summary(fun.y=mean,  geom='point', size=2.5) +
      stat_summary(fun.y=mean,  geom='line') +
      stat_summary(fun.data=mean_sdl, geom='linerange', alpha=0.4) +
      # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
      scale_color_manual("method",values=col) +
      scale_shape_manual("method",values=as.numeric(shape)) +
      scale_x_log10(breaks=n) +
      ylab("distance") + ggtitle(main) +
      theme_bw()+ theme(legend.title=element_text(size=14), legend.text=element_text(size=12))
}

pdf("plots/intramethod_distances.pdf", width=7, height=5)
for (i in 1:length(dmat_all))
  print(dist_plot(dmat_all[[i]], ddf, names(dmat_all)[i]))
dev.off()


n <- yaml::yaml.load_file('code/helpers/subsamples.yml')
dist_plot_cor <- function(Xcor, X, metadf, main="") {
   Xdf <- bind_rows(
            left_join(unpack(Xcor, dtypes_cor), metadf, by='dtype'),
            left_join(unpack(X, dtypes_all), metadf, by='dtype')
          )
   ggplot(aes(x=n, y=dist, col=factor(method, levels=allord), group=group,
              linetype=association, shape=factor(method, levels=allord)), data=Xdf) +
    stat_summary(fun.y=mean,  geom='point', data= subset(Xdf,association!="cor"), size=2.5) +
    stat_summary(fun.y=mean,  geom='line') +
    stat_summary(fun.data=mean_sdl, geom='linerange', alpha=0.4) +
   # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
    scale_color_manual("method",values=col) +
    scale_shape_manual("method",values=as.numeric(shape)) +
    scale_linetype_manual(values=c("dashed","solid","twodash"))+
    scale_x_log10(breaks=n) +
    ylab("distance") + ggtitle(main) +
    theme_bw()
}

pdf("plots/intramethod_distances_pearson_corshrink_propr.pdf", width=7, height=5)
for (i in 1:length(dmat_all_cor)) {
  nm <- names(dmat_all_cor)[i]
  print(dist_plot_cor(dmat_all_cor[[nm]], dmat_all[[nm]], ddf, nm))
}
dev.off()
