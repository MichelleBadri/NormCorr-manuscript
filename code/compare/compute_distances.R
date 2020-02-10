source('code/compare/metrics.R')

## Load in correlation estimates
est.list <- readRDS("data/RDS/est_all.RDS")
## Keep correlations/proportionalities and unlist 1 level
est.list <- unlist(lapply(est.list, function(x) x[grep("cor|rho", names(x))]), recursive=FALSE)

## Pairwise correlation matrix distances via parallelDist package
forst.all <- forstner.dist(est.list)
riemn.all <- riemannian.dist(est.list)
frob.all  <- frobenius.dist(est.list)
spect.all <- spectral.dist(est.list)
cmd.all   <- corrmat.dist(est.list)

## Save all as an RDS object
dmat_all <- list(forst=forst.all, riemn=riemn.all,
                 frob=frob.all,
                 spect=spect.all, cmd=cmd.all)
attr(dmat_all, 'datasets') <- names(est.list)

saveRDS(dmat_all, file="data/RDS/dmat_all.RDS")
