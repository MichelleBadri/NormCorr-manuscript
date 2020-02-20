suppressPackageStartupMessages(library(vsn))

# Read in normalized data matrices
est.list <- readRDS("data/RDS/est_all.RDS")
est.list <- est.list$`1_9000`
est.list_ag <-est.list[grep("amgut", names(est.list))]

est.list_ag$data.log2 <-t(log2(t(est.list_ag$amgut.raw+1)))

pdf("plots/meanSDplots.pdf", width=5.5,height=4)
for (i in 1:length(est.list_ag)) {
msd<-meanSdPlot(as.matrix(t(est.list_ag[[i]])),plot =FALSE)
print(msd$gg+geom_line(aes_string(x = "x", y = "y"),data = data.frame(x = msd[[1]], y = msd$sd), color="red", size=1.5) +theme_light() +ggtitle(names(est.list_ag)[i]))
}
dev.off()
