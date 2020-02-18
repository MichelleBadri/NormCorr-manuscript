suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggridges))
library(ggplot2)

source('code/correlation/plot_helpers.R')

## Read in full and shuffled correlation estimates
est.full  <- readRDS("data/RDS/est_all.RDS")
est.shuff <- readRDS("data/RDS/est_shuff.RDS")
## Compare 1 small and 1 large subsample
est.full <- est.full[names(est.shuff)]

## Keep (upper triangle of) correlation/proportionality estimates
unpackcorr <- function(x) lapply(x[grep("cor|rho", names(x))], SpiecEasi:::triu)
est.full  <- unlist(lapply(est.full, unpackcorr), recursive=FALSE)
est.shuff <- unlist(lapply(est.shuff, unpackcorr), recursive=FALSE)

fulldf  <- list2df(est.full)
shuffdf <- list2df(est.shuff)

shuffdf <- shuffdf[shuffdf$n == "50",]
fulldfsubsmall <- fulldf[fulldf$n=="50",]
fulldfsublarge <- fulldf[fulldf$n=="9000",]
shuffdf$group <- "n=50 shuffled" ## Titles for facets in plot
fulldfsubsmall$group<- "n=50"
fulldfsublarge$group<- "n=9000"

combdf <- rbind(shuffdf,fulldfsubsmall,fulldfsublarge)

pdf("plots/correlation_density.pdf", width=8, height=5)
plotDensity(combdf,levels=mainord)
dev.off()

## Change subset to plot supplementary ggridges plot
pdf("plots/correlation_density_supp.pdf", width=8, height=10)
plotDensity(combdf, levels=suppord)
dev.off()
