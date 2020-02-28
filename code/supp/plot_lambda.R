library(ggplot2)
library(reshape)
est.list <- readRDS("data/RDS/est_all.RDS")
## Keep 1 largest estimate

lamlist <- lapply(est.list, function(est) unlist(lapply(est, attr,"lambda")))

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
## order of methods for plotting legends
allord <- yaml::yaml.load_file('code/helpers/data_order.yml')[['allnorm']]


n <- yaml::yaml.load_file('code/helpers/subsamples.yml')

pdf("plots/lambda_shrinkage.pdf", width=7, height=4)
ggplot(aes(x=int, y=value, col=factor(method, levels=allord), shape=factor(method, levels=allord)), data=lamlistlong) +
  stat_summary(fun.y=mean,  geom='point', size=2.5,alpha=0.8) +
  stat_summary(fun.y=mean,  geom='line') +
  stat_summary(fun.data=mean_sdl, geom='linerange') +
  # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
  scale_color_manual("method",values=col) +
  scale_shape_manual("method",values=as.numeric(shape)) +
  scale_x_log10(breaks=n) +
  ylab("lambda") +
  theme_bw()
dev.off()


lambdavar <- lapply(est.list, function(est) unlist(lapply(est, attr,"lambda.var")))

rholamda<- data.frame(cbind(int=sapply(strsplit(names(est.list),split="\\_"), `[`, 2),lambda.var=lambdavar,lambda.cor=as.numeric(as.character(data.frame(lamlistwide)$rhoshrink))))
rholamda<- melt(data.frame(rholamda),id.vars = "int")
rholamda$int <- as.numeric(as.character(rholamda$int))
rholamda$value <- as.numeric(as.character(rholamda$value))

n <- yaml::yaml.load_file('code/helpers/subsamples.yml')

pdf("plots/lambda_rhoshrink.pdf", width=7, height=4)
ggplot(aes(x=int, y=value, col=variable, shape=variable), data=rholamda) +
  stat_summary(fun.y=mean,  geom='point', size=2.5) +
  stat_summary(fun.y=mean,  geom='line') +
  stat_summary(fun.data=mean_sdl, geom='linerange', alpha=0.4) +
  # stat_summary(fun.data=mean_sdl, geom='smooth',se=TRUE, alpha=0.15)+
  scale_color_manual("shrink type",values=c("#6F1E51" ,"#6F1E51" )) +
  scale_shape_manual("shrink type",values=c(1,16)) +
  scale_x_log10(breaks=n) +
  ylab("lambda") +
  theme_bw()
dev.off()
