suppressPackageStartupMessages(library(igraph))

igr_li <- readRDS("data/RDS/igraph_full.RDS")

## subset methods for final consensus
subset <- yaml::yaml.load_file('code/helpers/data_subsets.yml')[['consensus']]
igr_li <- igr_li[subset]

## extract signed edge list from igraph list
get_signed_edges <- function(igr_li) apply(cbind(get.edgelist(igr_li), E(igr_li)$sign), 1, paste, collapse="-")

edge_li <- lapply(igr_li, get_signed_edges)

## UpsetR plot
library(UpSetR)
pdf('plots/relevance_nets_upset.pdf', width=5, height=3)
out <- upset(fromList(edge_li), keep.order=TRUE,
            line.size=1)
out$Sizes[2] <- NULL
tryCatch(
print(out),
error=function(e) NULL
)
dev.off()

source('code/networks/plot_helpers.R')
## Create and plot venn diagram
pdf("plots/venndiagram.pdf", width=5.5, height=5)
do.call(draw.quad.venn, as.list(plot_venn_d(igr_li)))
dev.off()

## Create and plot consensus network
## Merge all graphs / attributes
igr_consens <- gr_intersect(igr_li)

pdf("plots/consensus_network.pdf", width=6, height=6)
plot_phy_net(igr_consens, "")
dev.off()
