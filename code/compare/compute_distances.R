source('code/compare/metrics.R')

## Load in correlation estimates
est.list <- readRDS("data/RDS/est_all.RDS")

## Calculate pearson correlation from amgut estimates
list.raw <- unlist(lapply(est.list, function(x) x[grep("amgut", names(x))]), recursive=FALSE)
cor.list <- lapply(list.raw,cor)

names(cor.list)<- paste(sapply(strsplit(names(cor.list),split="\\."), `[`, 1), sapply(strsplit(names(cor.list),split="\\."), `[`, 3),"cor",sep=".")

## Keep correlations/proportionalities and unlist 1 level
est.list <- unlist(lapply(est.list, function(x) x[grep("cor|rho", names(x))]), recursive=FALSE)

## Combine with corshrink and proportionality
est.list.cor <- c(est.list,cor.list)

#saveRDS(est.list.cor, file="data/RDS/est.list.corr.RDS")
#est.list.cor <- readRDS(file="data/RDS/est.list.corr.RDS")
## Pairwise correlation matrix distances via parallelDist package
frob.all  <- frobenius.dist(est.list)
spect.all <- spectral.dist(est.list)
cmd.all  <- corrmat.dist(est.list)
#forst.all <- forstner.dist(est.list)
#riemn.all <- riemannian.dist(est.list)

frob.all.cor  <- frobenius.dist(est.list.cor)
cmd.all.cor   <- corrmat.dist(est.list.cor)

## Save distances between corshrink/propr as an RDS object
dmat_all <- list(#forst=forst.all_cor, riemn=riemn.all_cor,
                     frob=frob.all,
                     spect=spect.all, 
                     cmd=cmd.all)
attr(dmat_all, 'datasets') <- names(est.list)

saveRDS(dmat_all, file="data/RDS/dmat_all.RDS")


## Save distances between pearson/corshrink/propr as an RDS object
dmat_all_cor <- list(#forst=forst.all_cor, riemn=riemn.all_cor,
                      frob=frob.all.cor,
                      cmd=cmd.all.cor)
attr(dmat_all_cor, 'datasets') <- names(est.list.cor)

saveRDS(dmat_all_cor, file="data/RDS/dmat_all_cor.RDS")
