## Minimal package requirements
library(SpiecEasi)
suppressPackageStartupMessages(library(phyloseq))
library(corpcor)
suppressPackageStartupMessages(library(dplyr))
library(Wrench)
suppressPackageStartupMessages(library(propr))

## Generate |sizes|*rep indices as a list for subsampling data
generate_indices <- function(n, sizes, reps) {
  imat <- expand.grid(rep=1:reps, size=sizes)
  list(imat, indices=lapply(imat$size, function(s) sample(n, size=s, replace=FALSE)))
}

## Apply norm/corr to subsample of OTU data
norm_corr <- function(ind, OTUDat) {
  ## Choose a random subsample
  OTU <- OTUDat[,ind]
  depths <- colSums(OTU)

  # Total sum scaling
  amgut_tss   <- t(OTU+1)/depths #Column-wise matrix normalization
  # CLR
  amgut_clr   <- t(SpiecEasi::clr(OTU+1, 2))
  # CSS
  amgut_css   <- suppressMessages(t(metagenomeSeq::cumNormMat(OTU+1, p=.5)))
  # RLE (DESeq)
  amgut_rle   <- t(DESeq.matrix(OTU+1))
  # VST (DESeq2)
  amgut_vst   <- t(VST(OTU+1))
  # COM
  mindepth  <- min(depths)
  amgut_com <- floor(t(OTU) * (mindepth/depths))
  # WREN
  amgut_wren <- suppressWarnings(wrench(OTU, condition=rep(2, dim(OTU)[2])))
  amgut_wren <- t(sweep(OTU,2,amgut_wren$nf,"/"))
  # RHO (Propr)
  suppressMessages(suppressWarnings(
  rho <- propr(t(OTU+1), metric = "rho", ivar='clr', symmetrize=TRUE, p=0)
  ))
  amgut_rhoprop <- (rho@matrix+t(rho@matrix))/2
  # RHO_SHRINK
  amgut_rhoshrink <- rho_shrink_est(OTU+1, verbose=FALSE)
  # RAW counts
  amgut_raw   <- t(OTU)
  # ASINH
  amgut_asinh <- t(scale(asinh(OTU+1), center=TRUE, scale=FALSE))

  ## Apply cor.shrink from corpcor to everything except rho and rhoshrink
  tss_cor  <- cor.shrink(amgut_tss, verbose=FALSE)
  clr_cor  <- cor.shrink(amgut_clr, verbose=FALSE)
  css_cor  <- cor.shrink(amgut_css, verbose=FALSE)
  rle_cor  <- cor.shrink(amgut_rle, verbose=FALSE)
  vst_cor  <- cor.shrink(amgut_vst, verbose=FALSE)
  com_cor  <- cor.shrink(amgut_com, verbose=FALSE)
  wren_cor <- cor.shrink(amgut_wren, verbose=FALSE)
  raw_cor  <- cor.shrink(amgut_raw, verbose=FALSE)
  asinh_cor <- cor.shrink(amgut_asinh, verbose=FALSE)

  ## convert all to class: matrix
  class(tss_cor)  <- "matrix"
  class(clr_cor)  <- "matrix"
  class(css_cor)  <- "matrix"
  class(rle_cor)  <- "matrix"
  class(vst_cor)  <- "matrix"
  class(com_cor)  <- "matrix"
  class(wren_cor) <- "matrix"
  class(raw_cor)  <- "matrix"
  class(amgut_rhoprop)  <- "matrix"
  class(amgut_rhoshrink) <- "matrix"
  class(asinh_cor)  <- "matrix"

  # save to an est object for saving for later (when testN > 2000, VST can take hours)
  est <- list(size= length(ind),
              ind = ind,
              raw.table=amgut_raw,
              amgut.tss= amgut_tss,
              amgut.clr= amgut_clr,
              amgut.css= amgut_css,
              amgut.rle= amgut_rle,
              amgut.vst= amgut_vst,
              amgut.com= amgut_com,
              amgut.wren= amgut_wren,
              amgut.raw= amgut_raw,
              amgut.asinh = amgut_asinh,
              tss.corshrink= tss_cor,
              clr.corshrink= clr_cor,
              css.corshrink= css_cor,
              rle.corshrink= rle_cor,
              vst.corshrink= vst_cor,
              com.corshrink= com_cor,
              wren.corshrink= wren_cor,
              raw.corshrink= raw_cor,
              asinh.corshrink= asinh_cor,
              rhoprop = amgut_rhoprop,
              rhoshrink = amgut_rhoshrink)
  return(est)
}
