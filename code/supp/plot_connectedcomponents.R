library(ggplot2)
library(reshape)
est.list <- readRDS("data/RDS/est_all.RDS")
## Keep 1 largest estimate
est <- est.list$'1_9000'
est.cor.rho <-est[grep("cor|rho", names(est))]
#rm(est.list)

## Load in data
suppressPackageStartupMessages(library(phyloseq))
ag.filt3 <- readRDS("data/RDS/ag.filt3.RDS")

maxeig=16

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

spectral.cluster.gap <- function(R, k1=3, k2='eigengap', kmax=20) {
  # R => Correlation/proportionality matrix
  # k1 => nearest neighbors
  # k2 => clustering strategy or k2 clusters
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
  return(rev(evL$values)[1:(kmax-1)])
}

cl_gap <- parallel::mclapply(est.cor.rho, spectral.cluster.gap,
                                k2="eigengap",kmax=maxeig)

gap <- cbind(eigen=seq(1:15),do.call(cbind, cl_gap))
eigs_melt<-melt.data.frame(data.frame(gap), id.vars = "eigen")

## Load method colors and shapes as a named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))
shape   <- unlist(yaml::yaml.load_file('code/helpers/pointshape.yml'))

eigs_melt$variable <- sapply(strsplit(as.character(eigs_melt$variable) ,split="\\."), `[`, 1)


pdf('plots/spectral_kcomponents.pdf')
lines <- ggplot(data=eigs_melt, aes(x=eigen, y=value,group=variable,colour=variable,shape=variable))+ geom_line(aes(y=value,colour=factor(variable)),alpha=0.6,size=1.3) + geom_point(size=2.6)+scale_color_manual(values = col) +xlab("Eigenvalue Rank") +ylab("Value")+scale_shape_manual(values=as.numeric(shape))
lines +theme_linedraw() +theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1))+ theme(legend.title=element_blank())  
dev.off()
