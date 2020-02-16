library(ggplot2)
library(reshape)
est.list <- readRDS("data/RDS/est_all.RDS")
## Keep 1 largest estimate

lamlist <- c()
for(i in 1:length(est.list)){
lamlist[[i]]<- unlist(lapply(est.list[[i]], attr,"lambda"))
}

lamlistwide<- do.call(rbind, lamlist)
lamlistwide<- cbind(int=sapply(strsplit(names(est.list),split="\\_"), `[`, 2),lamlistwide)
colnames(lamlistwide) <- sapply(strsplit(colnames(lamlistwide) ,split="\\."), `[`, 1)
lamlistlong<- melt(data.frame(lamlistwide),id.vars = "int")
lamlistlong$int <- as.numeric(as.character(lamlistlong$int))
lamlistlong$value <- as.numeric(as.character(lamlistlong$value))
lamlistlong <- lamlistlong[lamlistlong$variable!="rhoshrink",]
colnames(lamlistlong)[2] <- "method"
## method color named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))
## method shape named vector
shape   <- unlist(yaml::yaml.load_file('code/helpers/pointshape.yml'))
n <- c(25, 50, 100, 200, 350, 700, 1800, 4300, 9000)

pdf("plots/lambda_shrinkage.pdf", width=7, height=4)
ggplot(aes(x=int, y=value, col=method,shape=method), data=lamlistlong) +
  stat_summary(fun.y=mean,  geom='point', size=2.5) +
  stat_summary(fun.y=mean,  geom='line') +
  stat_summary(fun.data=mean_sdl, geom='linerange', alpha=0.4) +
  # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
  scale_color_manual(values=col) +
  scale_shape_manual(values=as.numeric(shape)) +
  scale_x_log10(breaks=n) +
  ylab("lambda") +
  theme_bw()
dev.off()


lambdavar <- c()
for(i in 1:length(est.list)){
  lambdavar[[i]]<- unlist(lapply(est.list[[i]], attr,"lambda.var"))
}


rholamda<- data.frame(cbind(int=sapply(strsplit(names(est.list),split="\\_"), `[`, 2),lambda.var=lambdavar,lambda.cor=as.numeric(as.character(data.frame(lamlistwide)$rhoshrink))))
rholamda<- melt(data.frame(rholamda),id.vars = "int")
rholamda$int <- as.numeric(as.character(rholamda$int))
rholamda$value <- as.numeric(as.character(rholamda$value))

pdf("plots/lambda_rhoshrink.pdf", width=7, height=4)
n <- c(25, 50, 100, 200, 350, 700, 1800, 4300, 9000)
ggplot(aes(x=int, y=value, col=variable,shape=variable), data=rholamda) +
  stat_summary(fun.y=mean,  geom='point', size=2.5) +
  stat_summary(fun.y=mean,  geom='line') +
  stat_summary(fun.data=mean_sdl, geom='linerange', alpha=0.4) +
  # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
  scale_color_manual(values=c("#6F1E51" ,"#6F1E51" )) +
  scale_shape_manual(values=c(1,16)) +
  scale_x_log10(breaks=n) +
  ylab("lambda") + ggtitle("covariance shrinkage for rho shrink") +
  theme_bw()
dev.off()