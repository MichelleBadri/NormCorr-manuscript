suppressPackageStartupMessages(library(dplyr))

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


## Load method colors as a named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))
library(ggplot2)
library(ggridges)

## Keep 4 methods for main ggridges plot
subset <- yaml::yaml.load_file('code/helpers/data_subsets.yml')[['correlation']]
subset <- stringr::str_split_fixed(subset, "_|\\.", 4)[,1]

plotDensity <- function(df, main="") {
  ggplot(aes(x=value, y=factor(method, levels=rev(subset)), fill=method), data=df %>% filter(method %in% subset)) +
      facet_wrap(~n) +
      geom_density_ridges(scale=1.5) +
      scale_fill_manual(values=col) +
      guides(fill=FALSE) + ylab("") +
      ## limit scales?
      scale_x_continuous(limits=c(-.5,.5)) +
      ggtitle(main) +
      theme_ridges(center_axis_labels = TRUE)
}


pdf("plots/correlation_density.pdf", width=6, height=4)
plotDensity(fulldf)
plotDensity(shuffdf, "Shuffled")
dev.off()
