## Write moments for density plots
est.full  <- readRDS("data/RDS/est_all.RDS")
est.shuff <- readRDS("data/RDS/est_shuff.RDS")


ind <- grep("cor|rho", names(est.full[[1]]))
est_50 <- est.full$`1_50`[ind]
est_9000 <- est.full$`1_9000`[ind]
est_50shuf <- est.shuff$`1_50`[ind]

rm(est.full)
rm(est.shuff)


for (i in 1:length(names(est_50))) {
#  out <- eval(parse(text=paste("est_50$", names(est_50)[i], sep="")))
  out <- est_50[[i]]
  out_utri <- out[upper.tri(out, diag=FALSE)]
#  out_utri <- out_utri[out_utri!=0]
  mean <- paste("mean =", round(mean(out_utri), digits=4))
  var <- paste("variance =", round(var(out_utri), digits=4))
  skew <- paste("skewness =", round(e1071::skewness(out_utri), digits=3))
  kurt <- paste("kurtosis =", round(e1071::kurtosis(out_utri), digits=3))
  cat(paste("\n",names(est_50)[i],mean,var,skew,kurt ,sep="\n"))
}


for(i in 1:length(names(est_9000))){
  out <- est_9000[[i]]
  out_utri <- out[upper.tri(out, diag=FALSE)]
  #  out_utri <- out_utri[out_utri!=0]
  mean <- paste("mean =", round(mean(out_utri), digits=4))
  var <- paste("variance =", round(var(out_utri), digits=4))
  skew <- paste("skewness =", round(e1071::skewness(out_utri), digits=3))
  kurt <- paste("kurtosis =", round(e1071::kurtosis(out_utri), digits=3))
  cat(paste("\n",names(est_9000)[i],mean,var,skew,kurt ,sep="\n"))
}


for(i in 1:length(names(est_50shuf))){
  out <- est_50shuf[[i]]
  out_utri <- out[upper.tri(out, diag=FALSE)]
  #  out_utri <- out_utri[out_utri!=0]
  mean <- paste("mean =", round(mean(out_utri), digits=4))
  var <- paste("variance =", round(var(out_utri), digits=4))
  skew <- paste("skewness =", round(e1071::skewness(out_utri), digits=3))
  kurt <- paste("kurtosis =", round(e1071::kurtosis(out_utri), digits=3))
  cat(paste("\n",names(est_50shuf)[i],mean,var,skew,kurt ,sep="\n"))
}
