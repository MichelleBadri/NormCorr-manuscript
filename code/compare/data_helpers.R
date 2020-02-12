library(stringr)
suppressPackageStartupMessages(library(dplyr))

## Read in distance matrices, method colors, metadata, etc...

dmat_all <- readRDS(file="data/RDS/dmat_all_cor.RDS")

## method color named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))

## method shape named vector
shape   <- unlist(yaml::yaml.load_file('code/helpers/pointshape.yml'))

## Dataset types
datasets <- attr(dmat_all, 'datasets')
dtypes <- stringr::str_split_fixed(datasets, "_", 2)[,2]
coldf <- data.frame(method=names(col), col=col, stringsAsFactors=FALSE)
shapedf <- data.frame(method=names(shape), shape=shape, stringsAsFactors=FALSE)
ddf <- stringr::str_split_fixed(datasets, "_|\\.", 4) %>%
          as.data.frame(stringsAsFactors=FALSE) %>%
          transmute(dtype=dtypes,
                 rep=as.numeric(V1),
                   n=as.numeric(V2),
                 method=V3,
                 association=V4,
                 group=paste(V3,V4,sep="")) %>%
          left_join(coldf, by = "method")
ddf<- left_join(shapedf,ddf, by = "method")
ddf[ddf$method=="rhoprop",]$association <- "propr"
ddf[ddf$method=="rhoshrink",]$association <- "corshrink"
rm(coldf)
rm(shapedf)
