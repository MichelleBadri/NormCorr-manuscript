## Write moments for density plots
est.full  <- readRDS("data/RDS/est_all.RDS")
est.shuff <- readRDS("data/RDS/est_shuff.RDS")

est_50 <- est.full$`1_50`[13:23]
est_9000 <- est.full$`1_9000`[13:23]
est_50shuf <- est.shuff$`1_50`[13:23]

rm(est.full)
rm(est.shuff)


for(i in 1:length(names(est_50))){
  out <- eval(parse(text=paste("est_50$", names(est_50)[i], sep="")))
  mean <- paste("mean =", round(mean(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=4))
  var <- paste("variance =", round(var(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=4))
  skew <- paste("skewness =", round(e1071::skewness(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=3))
  kurt <- paste("kurtosis =", round(e1071::kurtosis(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=3))
  cat(paste("\n",names(est_50)[i],mean,var,skew,kurt ,sep="\n"))
}


for(i in 1:length(names(est_9000))){
out <- eval(parse(text=paste("est_9000$", names(est_50)[i], sep="")))
mean <- paste("mean =", round(mean(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=4))
var <- paste("variance =", round(var(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=4))
skew <- paste("skewness =", round(e1071::skewness(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=3))
kurt <- paste("kurtosis =", round(e1071::kurtosis(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=3))
cat(paste("\n",names(est_50)[i],mean,var,skew,kurt ,sep="\n"))
}


for(i in 1:length(names(est_50shuf))){
  out <- eval(parse(text=paste("est_50shuf$", names(est_50)[i], sep="")))
  mean <- paste("mean =", round(mean(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=6))
  var <- paste("variance =", round(var(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=6))
  skew <- paste("skewness =", round(e1071::skewness(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=3))
  kurt <- paste("kurtosis =", round(e1071::kurtosis(out[upper.tri(out, diag=FALSE)][out[upper.tri(out, diag=FALSE)]!=0]), digits=3))
  cat(paste("\n",names(est_50)[i],mean,var,skew,kurt ,sep="\n"))
}
