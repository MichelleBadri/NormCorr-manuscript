source("code/helpers/norm_functions.R")
source("code/correlation/draw_samples.R")

ag.filt3 <- readRDS("data/RDS/ag.filt3.RDS")

set.seed(10010)
#OTUDat <- ag.filt3@otu_table@.Data
## Permute counts between samples (columns)
OTUDat <- t(apply(ag.filt3@otu_table@.Data, 1, sample))
colnames(OTUDat) <- sample_names(ag.filt3)

## 1 small and 1 large subsample size and a single replicate
testN <- c(50, 9000)
reps  <- 1

## Generate subsamples using a fixed seed
set.seed(1)
tmp <- generate_indices(ncol(OTUDat), testN, reps)
sizemat <- tmp[[1]]
indices <- tmp[[2]]

## Use random indices to normalize counts and estimate correlations
est.list <- parallel::mclapply(indices, norm_corr, OTUDat=OTUDat,
                        mc.cores=parallel::detectCores(),
                        mc.preschedule=FALSE)
## List name is a combination of parameter matrix
names(est.list) <- apply(sizemat, 1, paste, collapse="_")

saveRDS(est.list, file="data/RDS/est_shuff.RDS")

