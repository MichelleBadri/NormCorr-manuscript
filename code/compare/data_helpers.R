library(stringr)
suppressPackageStartupMessages(library(dplyr))

## Read in distance matrices, method colors, metadata, etc...


## Load distances with pearson, corshrink and propr
dmat_all_cor <- readRDS(file="data/RDS/dmat_all_cor.RDS")

## method color named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))

## method shape named vector
shape   <- unlist(yaml::yaml.load_file('code/helpers/pointshape.yml'))

## Dataset types
datasets_cor <- attr(dmat_all_cor, 'datasets')
dtypes_cor <- stringr::str_split_fixed(datasets_cor, "_", 2)[,2]
coldf <- data.frame(method=names(col), col=col, stringsAsFactors=FALSE)
shapedf <- data.frame(method=names(shape), shape=shape, stringsAsFactors=FALSE)
ddf_cor <- stringr::str_split_fixed(datasets_cor, "_|\\.", 4) %>%
          as.data.frame(stringsAsFactors=FALSE) %>%
          transmute(dtype=dtypes_cor,
                 rep=as.numeric(V1),
                 n=as.numeric(V2),
                 method=V3,
                 association=V4,
                 group=paste(V3,V4,sep="")) %>%
          left_join(coldf, by = "method")
ddf_cor<- left_join(shapedf,ddf_cor, by = "method")
ddf_cor[ddf_cor$method=="rhoprop",]$association <- "propr"
ddf_cor[ddf_cor$method=="rhoshrink",]$association <- "corshrink"



## Load distances with corshrink and propr
dmat_all <- readRDS(file="data/RDS/dmat_all.RDS")

## Dataset types
datasets <- attr(dmat_all, 'datasets')
dtypes <- stringr::str_split_fixed(datasets, "_", 2)[,2]
coldf <- data.frame(method=names(col), col=col, stringsAsFactors=FALSE)
ddf <- stringr::str_split_fixed(datasets, "_|\\.", 4) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  transmute(dtype=dtypes,
            rep=as.numeric(V1),
            n=as.numeric(V2),
            method=V3,
            association=V4) %>%
  left_join(coldf, by = "method")
ddf<- left_join(shapedf,ddf, by = "method")
ddf[ddf$method=="rhoprop",]$association <- "propr"
ddf[ddf$method=="rhoshrink",]$association <- "corshrink"

rm(coldf)
rm(shapedf)



