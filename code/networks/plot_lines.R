suppressPackageStartupMessages(library(igraph))
library(stringr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


igr_li_allsub <- readRDS("data/RDS/igraph_subsets.RDS")

## extract community analysis measures from graph objects
Assortativity_Genus <- unlist(lapply(igr_li_allsub, get.graph.attribute,"assortativity_genus"))
Maximum_Modularity <- unlist(lapply(igr_li_allsub, get.graph.attribute,"max_modularity"))

datasets <- attr(igr_li_allsub, 'names')
assort_ddf <- data.frame(datasets,Assortativity_Genus)
assort_ddf$datasets <- as.character(assort_ddf$datasets)

maxmod_ddf <- data.frame(datasets,Maximum_Modularity)
maxmod_ddf$datasets <- as.character(maxmod_ddf$datasets)


## method color named vector
col <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))

## method shape named vector
shape <- unlist(yaml::yaml.load_file('code/helpers/pointshape.yml'))

## order of methods for plotting legends
allord <- yaml::yaml.load_file('code/helpers/data_order.yml')[['allnorm']]

## Set up long data for ggplot2 plotting
dtypes <- stringr::str_split_fixed(datasets, "_", 2)[,2]
coldf <- data.frame(method=names(col), col=col, stringsAsFactors=FALSE)
shapedf <- data.frame(method=names(shape), shape=shape, stringsAsFactors=FALSE)
ddf_cor <- stringr::str_split_fixed(datasets, "_|\\.", 4) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  transmute(datasets=datasets,
            dtype=dtypes,
            rep=as.numeric(V1),
            n=as.numeric(V2),
            method=V3,
            association=V4,
            group=paste(V3,V4,sep="")) %>%
  left_join(coldf, by = "method")
ddf_cor<- left_join(shapedf,ddf_cor, by = "method")
assort_ddf<- left_join(assort_ddf,ddf_cor, by = "datasets")
maxmod_ddf<- left_join(maxmod_ddf,ddf_cor, by = "datasets")
rm(coldf)
rm(shapedf)


n <- yaml::yaml.load_file('code/helpers/subsamples.yml')
pdf("plots/relevance_nets_communityanalysis_lineplots.pdf", width=7, height=4)
Xdf <-assort_ddf
  ggplot(aes(x=n, y=Assortativity_Genus, col=factor(method, levels=allord),
            shape=factor(method, levels=allord)), data=Xdf) +
    stat_summary(fun.y=mean,  geom='point', size=2.5,alpha=0.8) +
    stat_summary(fun.y=mean,  geom='line') +
    stat_summary(fun.data=mean_sdl, geom='linerange') +
    # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
    scale_color_manual("method",values=col) +
    scale_shape_manual("method",values=as.numeric(shape)) +
    scale_x_log10(breaks=n) +
    theme_bw()

Xdf <- maxmod_ddf
  ggplot(aes(x=n, y=Maximum_Modularity, col=factor(method, levels=allord),
            shape=factor(method, levels=allord)), data=Xdf) +
    stat_summary(fun.y=mean,  geom='point', size=2.5,alpha=0.8) +
    stat_summary(fun.y=mean,  geom='line') +
    stat_summary(fun.data=mean_sdl, geom='linerange') +
    # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
    scale_color_manual("method",values=col) +
    scale_shape_manual("method",values=as.numeric(shape)) +
    scale_y_continuous(limits=c(0.4,1)) +
    scale_x_log10(breaks=n) +
    theme_bw()
dev.off()
