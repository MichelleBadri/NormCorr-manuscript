suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(circlize))

est.list <- readRDS("data/RDS/est_all.RDS")
## Keep 1 large estimate
est <- est.list$'1_9000'
rm(est.list)

hier.cluster <- function(R) {
  D    <- sqrt((1-R)/2)
  Dmat <- as.dist(D)
  hclust(Dmat, method="ward.D2") #Ward D2 squares the distances
}


subset <- yaml::yaml.load_file('code/helpers/data_subsets.yml')[['correlation']]
hcli <- parallel::mclapply(est[subset], hier.cluster,
                         mc.cores=parallel::detectCores(),
                         mc.preschedule=FALSE)

source('code/clustering/plot_helpers.R')

pdf('plots/hierarchical_clustering.pdf')
for (i in 1:length(hcli)) {
  plot_circos(hcli[[i]], names(hcli)[i])
}
dev.off()

## Plot supplement
supp <- yaml::yaml.load_file('code/helpers/data_subsets.yml')[['supplement']]
 
hcli_supp <- parallel::mclapply(est[supp], hier.cluster,
                           mc.cores=parallel::detectCores(),
                           mc.preschedule=FALSE)

source('code/clustering/plot_helpers.R')

pdf('plots/hierarchical_clustering_supp.pdf')
for (i in 1:length(hcli_supp)) {
  plot_circos(hcli_supp[[i]], names(hcli_supp)[i])
}
dev.off()
