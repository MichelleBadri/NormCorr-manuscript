suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

## load dmat_all, dataset df (ddf), color vector (col)
source("code/compare/data_helpers.R")

## compute all embeddings via isometric MDS
isomds <- function(x) {set.seed(10010) ; MASS::isoMDS(x, k=3)$points}
mds_all <- parallel::mclapply(dmat_all, isomds,
                          mc.cores=parallel::detectCores(),
                          mc.preschedule=FALSE)
names(mds_all) <- names(dmat_all)

## Rotate embeddings w.r.t. the Frobenius dist embedding
sind <- !(names(mds_all)=="frob")
proc <- function(x) vegan::procrustes(mds_all$frob, x)$Yrot
mds_all[sind] <- lapply(mds_all[sind], proc)

## 3d plots
plot3d <- function(X, metadf, col) {
  stopifnot(nrow(X)==nrow(metadf))
  Xdf <- as.data.frame(X)
  plotly::plot_ly(x=Xdf[,1], y=Xdf[,2], z=Xdf[,3], color=metadf$method,
                  colors=col, sizes = c(50,200), size=metadf$n)
}


plot2d <- function(X, metadf, main="") {
  ## every other sample index
  n <- c(25, 100, 350, 1800, 9000)
  stopifnot(nrow(X)==nrow(metadf))
  ## Embeddings to dataframe for ggplot2
  Xdf <- cbind(as.data.frame(X), ddf)
  ## Plot embedding group centers for each method/sample size
  Xmean <- Xdf %>% group_by(method, n) %>% summarize(V1=mean(V1), V2=mean(V2))
  ggplot(aes(x=V1, y=V2, col=method, size=n), data=Xdf) +
    geom_point(alpha=.2) +
    scale_color_manual(values=col) +
    scale_size_continuous(breaks=n) +
    geom_point(data=Xmean, alpha=.8) +
    xlab("MDS1") + ylab("MDS2") + ggtitle(main) +
    theme_bw()
}

## Three options for projection matrices
## (estimated by hand-rotating plot in rgl::plot3d and calling rgl::par3d()$userMatrix)
A <- matrix(c(1,0,0,0,.8,-.6,0,.6,.8),nrow=3)
B <- matrix(c(.6,.1,.8,.8,-.04,-.6,.04,1,.1),nrow=3)
C <- matrix(c(.7,-.5,4,.7,.5,-.5,-.025,.7,.7), nrow=3)

pdf("plots/embeddings_2d.pdf")
for (i in 1:length(mds_all))
 print(plot2d(mds_all[[i]] %*% t(C), ddf, names(mds_all)[i]))
dev.off()
