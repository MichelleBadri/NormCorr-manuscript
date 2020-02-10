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
spectral.cluster <- function(R, k1=3, k2='eigengap', kmax=20) {
  # R => Correlation/proportionality matrix
  # k1 => nearest neighbors
  # k2 => clustering strategy or k2 clusters
  # kmax => maximum number of 'clusters' to consider for the eigengap method

   ## construct similarity matrix
  S <- 1-(sqrt(1-R)/2)
  A <- make.affinity(S , 2, 'or')
  ## degree matrix
  d <- apply(A, 1, sum)
  D <- diag(d)
  D_isqrt <- diag(1/sqrt(d))
  ## normalized, symmetric Laplacian
  L <- D_isqrt %*% (D - A) %*% D_isqrt
  evL <- eigen(L, symmetric=TRUE)
  if (k2 == "eigengap") {
    ## check for the largest gap in the first 20 eigenvalues
    k2 <- which.max((diff(rev(evL$values)[1:(kmax-1)])))
  } else if (k2 == "components") {
    zeroEvals <- which(rev(abs(evL$values)) < 1e-9)
    tmp <- rle(diff(zeroEvals))
    k2 <- tmp$lengths[1]+1
  }  else if (is.numeric(k2)) {
    ## k2 is the number of clusters
  }
  else {
    .NotYetImplemented("No other methods available")
  }

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

## Try 2 cluster detection methods
cl_egap20 <- parallel::mclapply(est[subset], spectral.cluster,
                         k2="eigengap",
                         mc.cores=parallel::detectCores(), mc.preschedule=FALSE)

cl_comp <- parallel::mclapply(est[subset], spectral.cluster,
                        k2="components",
                        mc.cores=parallel::detectCores(), mc.preschedule=FALSE)



source('code/clustering/plot_helpers.R')

pdf('plots/spectral_clustering_eigengap.pdf')
for (i in 1:length(cl_egap20)) {
  mean_purity <- mean(cluster_purity(cl_egap20[[i]]))
  title <- paste(names(cl_egap20)[i], signif(mean_purity,2), sep=": ")
  print(plot_kmeans(cl_egap20[[i]], title))
}
dev.off()

pdf('plots/spectral_clustering_components.pdf')
for (i in 1:length(cl_comp)) {
  mean_purity <- mean(cluster_purity(cl_comp[[i]]))
  title <- paste(names(cl_comp)[i], signif(mean_purity,2), sep=": ")
  print(plot_kmeans(cl_comp[[i]], title))
}
dev.off()
