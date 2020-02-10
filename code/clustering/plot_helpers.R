library(ggplot2)

ag.filt3 <- readRDS("data/RDS/ag.filt3.RDS")

fcol <- as.vector(yaml::yaml.load_file('code/helpers/fam_colors.yml'))
plot_kmeans <- function(km, main="") {
  dat <- data.frame(cluster=as.factor(km$cluster),
                    Rank5=ag.filt3@tax_table@.Data[,5])
  counttab <- prop.table(table(dat$cluster, dat$Rank5),1)

  ord <- names(sort(-colSums(counttab)))
  dat <- reshape2::melt(counttab)
  order2 <- rownames(counttab[order(-counttab[,"f__Ruminococcaceae"]),, drop=FALSE])
  dat <- dat[order(match(dat$Var2,ord)),]
  ggplot(dat, aes(x=factor(dat$Var1,levels=order2),
                  y=value,
                  fill=Var2,
                  label=Var2,
                  group=factor(dat$Var2,levels=rev(ord)))) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x="Cluster Size", y="Proportion") +
    scale_fill_manual(values = fcol) +
    coord_flip() + guides(fill=FALSE) +
    scale_x_discrete(labels=as.character(km$size[as.numeric(order2)]),
                     expand = c(0, 0)) +
    scale_y_continuous(limits=c(0,1.008), expand = c(0, 0)) +
    ggtitle(main) +
    theme_classic() +
    theme(axis.text.y = element_text(face="bold", size=12),
          axis.ticks.y=element_blank(),
          plot.title=element_text(face="bold",size=18))
}


rank4 <- ag.filt3@tax_table@.Data[,"Rank4"]
ocol <- unlist(yaml::yaml.load_file('code/helpers/order_colors.yml'))
plot_circos <- function(hc, main) {
  dend <- ladderize(as.dendrogram(hc))
  ocols <- ocol[rank4[order.dendrogram(dend)]]
  n <- length(hc$order)
  circos.par(cell.padding = c(0, 0, 0, 0))
  # only one sector
  circos.initialize(factors = "a", xlim = c(0, n))
  # maximum height of the trees
  max_height <- attr(dend, "height")
  k <- 1
  pfun1 <- function(x, y) {
    for (i in seq_len(n)) circos.points(x=i,y=0.0002,pch=20,col=ocols[i], cex=1.8)
  }
  pfun2 <- function(x, y) circos.dendrogram(dend, max_height = max_height)
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.5, panel.fun = pfun1)
  circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, panel.fun = pfun2)
  title(main)
}
