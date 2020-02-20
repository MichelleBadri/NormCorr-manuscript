est.list <- readRDS("data/RDS/est_all.RDS")
## Keep 1 largest estimate
est <- est.list$'1_9000'
rm(est.list)

## Load in data
suppressPackageStartupMessages(library(phyloseq))
ag.filt3 <- readRDS("data/RDS/ag.filt3.RDS")

## Create affinity matrix from nearest neighbor similarity matrix
make.affinity <- function(S, k=3, sym="or") {
  S <- S-diag(diag(S))
  if (k>=ncol(S)) {
    return(S)
  } else {
    thr <- apply(S, 1, function(s) s[order(s, decreasing=TRUE)[k]])
    A <- S>=thr
    if (sym == "or") A <- sign(A + t(A))
    else if (sym == "and") A <- sign(A * t(A))
    return(S*A)
  }
}

## Normalized spectral clustering, Ng, Jordan & Weiss (2002)
spectral.cluster <- function(R, k1=3, kmax=20) {
  # R => Correlation/proportionality matrix
  # k1 => nearest neighbors
  # kmax => maximum number of 'clusters' to consider for the eigengap method
  
  ## construct similarity matrix
  S <- 1-(sqrt((1-R)/2))
  A <- make.affinity(S , 2, 'or')
  ## degree matrix
  d <- apply(A, 1, sum)
  D <- diag(d)
  D_isqrt <- diag(1/sqrt(d))
  ## normalized, symmetric Laplacian
  L <- D_isqrt %*% (D - A) %*% D_isqrt
  evL <- eigen(L, symmetric=TRUE)
  zeroEvals <- which(rev(abs(evL$values)) < 1e-9)
  tmp <- rle(diff(zeroEvals))
  k2 <- tmp$lengths[1]+1
  T  <- evL$vectors[,1:k2]
  T  <- T/sqrt(rowSums(T^2))
  # define 0/0 := 0
  T[!is.finite(T)] <- 0
  set.seed(1000)
  km <- kmeans(T, centers=k2, nstart=5000)
  ## Return eigendecomposition for spectrum plots
  km$Leig <- evL
  ## Return Family counts per cluster
  km$counttab <- table(km$cluster, ag.filt3@tax_table@.Data[,5])
  return(km)
}


cluster_purity <- function(km) {
  ## Cluster 'purity' is neff/species
  ## Single family cluster == 1
  neff <- exp(vegan::diversity(km$counttab, "shannon"))
  maxneff <- vegan::specnumber(km$counttab)
  neff/maxneff
}

## subset methods for main plots
subset <- yaml::yaml.load_file('code/helpers/data_subsets.yml')[['correlation']]

cl_comp <- parallel::mclapply(est[subset], spectral.cluster,
                              mc.cores=parallel::detectCores(), mc.preschedule=FALSE)

names(cl_comp) <- sapply(strsplit(as.character(names(cl_comp)) ,split="\\."), `[`, 1)

source('code/clustering/plot_helpers.R')

pdf('plots/spectral_clustering_components.pdf', width=5.5, height=5)
for (i in 1:length(cl_comp)) {
  mean_purity <- mean(cluster_purity(cl_comp[[i]]))
  title <- paste(names(cl_comp)[i], signif(mean_purity,2), sep=": ")
  print(plot_kmeans(cl_comp[[i]], title))
}
dev.off()


supp <- yaml::yaml.load_file('code/helpers/data_subsets.yml')[['supplement']]

cl_comp_supp <- parallel::mclapply(est[supp], spectral.cluster,
                              mc.cores=parallel::detectCores(), mc.preschedule=FALSE)

names(cl_comp_supp) <- sapply(strsplit(as.character(names(cl_comp_supp)) ,split="\\."), `[`, 1)

pdf('plots/spectral_clustering_components_supp.pdf', width=5.5, height=5)
for (i in 1:length(cl_comp_supp)) {
  mean_purity <- mean(cluster_purity(cl_comp_supp[[i]]))
  title <- paste(names(cl_comp_supp)[i], signif(mean_purity,2), sep=": ")
  print(plot_kmeans(cl_comp_supp[[i]], title))
}
dev.off()