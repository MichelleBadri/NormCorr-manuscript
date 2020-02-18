suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggridges))


## Read in full and shuffled correlation estimates
est.full  <- readRDS("data/RDS/est_all.RDS")
est.shuff <- readRDS("data/RDS/est_shuff.RDS")
## Compare 1 small and 1 large subsample
est.full <- est.full[names(est.shuff)]

## Keep (upper triangle of) correlation/proportionality estimates
unpackcorr <- function(x) lapply(x[grep("cor|rho", names(x))], SpiecEasi:::triu)
est.full  <- unlist(lapply(est.full, unpackcorr), recursive=FALSE)
est.shuff <- unlist(lapply(est.shuff, unpackcorr), recursive=FALSE)

## Unlist correlations/ convert to data frames for ggplot2/ggridges
list2df <- function(est.list) {
  reshape2::melt(est.list) %>%
          rename(dataset=L1) %>%
          bind_cols(
            stringr::str_split_fixed(.$dataset, "_|\\.", 4) %>%
                      as.data.frame(stringsAsFactors=FALSE) %>%
                      mutate(rep=as.numeric(V1),
                               n=as.numeric(V2),
                             method=V3) %>%
                      select(rep, n, method)
          )
}
fulldf  <- list2df(est.full)
shuffdf <- list2df(est.shuff)
shuffdf <- shuffdf[shuffdf$n == "50",]
fulldfsubsmall <- fulldf[fulldf$n=="50",]
fulldfsublarge <- fulldf[fulldf$n=="9000",]
shuffdf$group <- "n=50 shuffled"
fulldfsubsmall$group<- "n=50"
fulldfsublarge$group<- "n=9000"

combdf <- rbind(shuffdf,fulldfsubsmall,fulldfsublarge)

## Load method colors as a named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))
library(ggplot2)
library(ggridges)

## Keep 4 methods for main ggridges plot
subset <- yaml::yaml.load_file('code/helpers/data_subsets.yml')[['correlation']]
subset <- stringr::str_split_fixed(subset, "_|\\.", 4)[,1]

plotDensity <- function(df, main="", levels=order) {
  ggplot(aes(x=value, y=factor(method, levels=rev(levels)), fill=method), data=df %>% filter(method %in% subset)) +
      facet_wrap(~group) + #,scales = "free_x") +
      geom_density_ridges(scale=1.3) +
      scale_fill_manual(values=col) +
      guides(fill=FALSE) + ylab("") +
      ## limit scales?
      scale_x_continuous(limits=c(-0.25,0.25), breaks=c(-0.2,-0.1,0,0.1,0.2)) + # comment out for free scales- see full tails
      ggtitle(main) +
      theme_ridges(font_size = 14,center_axis_labels = TRUE) + theme_bw()
}

mainord <- yaml::yaml.load_file('code/helpers/data_order.yml')[['main']]
suppord <- yaml::yaml.load_file('code/helpers/data_order.yml')[['supp']]

pdf("plots/correlation_density.pdf", width=6, height=4)
plotDensity(combdf,levels=mainord)
dev.off()

## Change subset to plot supplementary ggridges plot
subset<- unique(combdf[!combdf$method %in% subset,]$method)

pdf("plots/correlation_density_supp.pdf", width=8, height=10)
plotDensity(combdf, levels=suppord)
dev.off()
