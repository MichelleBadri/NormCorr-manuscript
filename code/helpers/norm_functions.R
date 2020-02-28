## DESeq.matrix function used to calculate RLE (see NormCorr RLE methods)
DESeq.matrix <- function(mat, c) {
  # compute geometric mean along row
  matt  <- mat
  matt[matt == 0] <- NA
  k_ref <- apply(matt, 1, function(x) exp(mean(log(x), na.rm=TRUE)))
  krefmat <- matrix(rep(k_ref,ncol(mat)), nrow=nrow(mat))
  s_hat   <- apply(matt/krefmat, 2, median, na.rm=TRUE)
  if (missing(c)) {
    fn <- function(c, s_hat) abs(sum(log(c*s_hat)))
    c  <- optimize(fn, interval=0:10, s_hat=s_hat)$minimum
  }
  s <- c * s_hat
  smat <- matrix(rep(s,nrow(mat)), ncol=ncol(mat), byrow=TRUE)
  mat/smat
}

library(corpcor)
library(SpiecEasi)
rho_shrink_est <- function(OTU, ...) {
  OTU_clr <- t(SpiecEasi::clr(OTU, 2))
  shrunk_cov<- cov.shrink(OTU_clr, ...)
  p <- ncol(OTU_clr)
  J <- matrix(rep(diag(shrunk_cov), p), p)
  rho <- 2*shrunk_cov/(J+t(J))
  (rho + t(rho)) / 2
}

suppressPackageStartupMessages(library(DESeq2))
VST <- function(OTU) {
  obj <- suppressMessages(
            DESeqDataSetFromMatrix(OTU,
              DataFrame(row.names = colnames(OTU)), ~1)) %>%
          estimateSizeFactors
  design(obj) <- ~1

  obj %>%
      estimateDispersionsGeneEst(quiet=TRUE, dispTol=1e-1) %>%
      estimateDispersionsFit(quiet=TRUE, "local") %>%
      getVarianceStabilizedData
}