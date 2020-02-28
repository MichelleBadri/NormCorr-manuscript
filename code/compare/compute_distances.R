source('code/compare/metrics.R')

## Load in correlation estimates
est.list <- readRDS("data/RDS/est_all.RDS")

## Calculate pearson correlation from amgut estimates
list.raw <- unlist(lapply(est.list, function(x) x[grep("amgut", names(x))]), recursive=FALSE)

## Keep correlations/proportionalities and unlist 1 level
est.list <- unlist(lapply(est.list, function(x) x[grep("cor|rho", names(x))]), recursive=FALSE)

# ## Combine with corshrink and proportionality
# est.list.cor <- c(est.list,cor.list)

#saveRDS(est.list.cor, file="data/RDS/est.list.corr.RDS")
#est.list.cor <- readRDS(file="data/RDS/est.list.corr.RDS")
## Pairwise correlation matrix distances via parallelDist package
frob.all  <- frobenius.dist(est.list)
spect.all <- spectral.dist(est.list)
cmd.all  <- corrmat.dist(est.list)
#forst.all <- forstner.dist(est.list)
#riemn.all <- riemannian.dist(est.list)

## Save distances between corshrink/propr as an RDS object
dmat_all <- list(frob=frob.all,
                 spect=spect.all,
                 cmd=cmd.all)
attr(dmat_all, 'datasets') <- names(est.list)

saveRDS(dmat_all, file="data/RDS/dmat_all.RDS")

## Compute distances on 'raw' correlation matrices
corNA <- function(x) {
##  convert NAs to 0s in case there are any zero stddev vars
  corX <- suppressWarnings(cor(x))
  corX[is.na(corX)] <- 0
  corX
}
cor.list <- parallel::mclapply(list.raw, corNA,
                    mc.cores=parallel::detectCores(), mc.preschedule=FALSE)
names(cor.list) <- paste0(gsub("amgut.", "", names(cor.list)), ".cor")

## We only need to calculate intra-method distances for raw correlations
cor.list.t <- split(cor.list,  sapply(strsplit(names(cor.list), "\\."), '[', 2))
frob.all.cor  <- lapply(cor.list.t, frobenius.dist)
cmd.all.cor   <- lapply(cor.list.t, corrmat.dist)

## Save distances between pearson/corshrink/propr as an RDS object
dmat_all_cor <- list(frob=frob.all.cor,
                      cmd=cmd.all.cor)
attr(dmat_all_cor, 'datasets') <- lapply(cor.list.t, names)

saveRDS(dmat_all_cor, file="data/RDS/dmat_all_cor.RDS")
