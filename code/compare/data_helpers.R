library(stringr)
suppressPackageStartupMessages(library(dplyr))

## Read in distance matrices, method colors, metadata, etc...

dmat_all <- readRDS(file="data/RDS/dmat_all.RDS")

## method color named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))

## Dataset types
datasets <- attr(dmat_all, 'datasets')
dtypes <- stringr::str_split_fixed(datasets, "_", 2)[,2]
coldf <- data.frame(method=names(col), col=col, stringsAsFactors=FALSE)
ddf <- stringr::str_split_fixed(datasets, "_|\\.", 4) %>%
          as.data.frame(stringsAsFactors=FALSE) %>%
          transmute(dtype=dtypes,
                 rep=as.numeric(V1),
                   n=as.numeric(V2),
                 method=V3) %>%
          left_join(coldf, by = "method")
rm(coldf)