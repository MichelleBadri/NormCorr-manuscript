source("code/helpers/norm_functions.R")
source("code/correlation/draw_samples.R")

ag.filt3 <- readRDS("data/RDS/ag.filt3.RDS")

OTUDat <- ag.filt3@otu_table@.Data

## Subsample sizes and replicates
#testN <- c(25,50,100,250,500,1000,1500,2000)
testN <- c(25, 50, 100, 200, 350, 700, 1800, 4300, 9000)
reps  <- 5

## Generate subsampling indexes using a fixed seed
set.seed(1)
tmp <- generate_indices(ncol(OTUDat), testN, reps)
sizemat <- tmp[[1]]
indices <- tmp[[2]]

## Use random indices to normalize counts and estimate correlations
est.list <- parallel::mclapply(indices, norm_corr, OTUDat=OTUDat,
                        mc.cores=parallel::detectCores(), mc.preschedule=FALSE)
## List name is a combination of parameter matrix
names(est.list) <- apply(sizemat, 1, paste, collapse="_")

saveRDS(est.list, file="data/RDS/est_all.RDS")
