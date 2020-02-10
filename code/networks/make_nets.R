suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(phyloseq))

ag.filt3 <- readRDS("data/RDS/ag.filt3.RDS")
est.list <- readRDS("data/RDS/est_all.RDS")

## use 1 large sample dataset
est <- est.list$'1_9000'
## keep only correlation/proportionality estimates
est <- est[grep("cor|rho", names(est))]
rm(est.list)

cor_to_graph <- function(R, nedges=2000, rseed=10010) {
  ## R => correlation matrix
  ## nedges => rank edges to keep
  ## rseed => random seed for layout

  ## convert abs correlations to edge rank matrix sans diagonal/ keep `nedges` edges
  rankR <- matrix(rank(-abs(R-diag(diag(R)))), dim(R))
  adj <- (rankR<=nedges*2)*R
  G <- graph.adjacency(adj, weighted=TRUE, mode='undirected')

  ## Weights are correlation magnitudes / keep sign data as edge attributes
  E(G)$sign   <- sign(E(G)$weight)
  E(G)$weight <- abs(E(G)$weight)

  ## Store greedy cluster membership as vertex attributes
  fgreedy <- fastgreedy.community(G, merges=TRUE, modularity=TRUE)
  E(G)$cross  <- crossing(fgreedy, G)
  V(G)$cmemb  <- fgreedy$membership

  ## Store Phylum and Family as vertex attributes
  V(G)$Rank2 <- ag.filt3@tax_table@.Data[V(G)$name,"Rank2"]
  V(G)$Rank5 <- ag.filt3@tax_table@.Data[V(G)$name,"Rank5"]

  ## Store FR layout as vertex attributes
  set.seed(rseed)
  l <- layout.fruchterman.reingold(G)
  V(G)$l1 <- l[,1]
  V(G)$l2 <- l[,2]
  G
}

igr_li <- parallel::mclapply(est, cor_to_graph,
                         mc.cores=parallel::detectCores(),
                         mc.preschedule=FALSE)

saveRDS(igr_li, file="data/RDS/igraph_full.RDS")
