library(stringr)
suppressPackageStartupMessages(library(dplyr))

# Read in distance matrices, method colors, metadata, etc...

## Load distances with pearson, corshrink and propr
dmat_all_cor <- readRDS(file="data/RDS/dmat_all_cor.RDS")
dmat_all <- readRDS(file="data/RDS/dmat_all.RDS")

## method color named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))

## method shape named vector
shape   <- unlist(yaml::yaml.load_file('code/helpers/pointshape.yml'))

## order of methods for plotting legends
allord <- yaml::yaml.load_file('code/helpers/data_order.yml')[['allnorm']]
#
## Dataset types
dtypes_all <- stringr::str_split_fixed(unique(attr(dmat_all, 'datasets')), "_", 2) [,2]
dtypes_cor <-  lapply(attr(dmat_all_cor, 'datasets'), function(x) stringr::str_split_fixed(unique(x), "_", 2)[,2])
dtypes <- c(dtypes_all, unlist(dtypes_cor))
coldf <- data.frame(method=names(col), col=col, stringsAsFactors=FALSE)
shapedf <- data.frame(method=names(shape), shape=shape, stringsAsFactors=FALSE)

datasets <- unique(c(attr(dmat_all, 'datasets'),
                   unlist(attr(dmat_all_cor, 'datasets'))))
ddf <- stringr::str_split_fixed(datasets, "_|\\.", 4) %>%
          as.data.frame(stringsAsFactors=FALSE) %>%
          transmute(dtype=dtypes,
                 rep=as.numeric(V1),
                 n=as.numeric(V2),
                 method=V3,
                 association=V4,
                 group=paste(V3,V4,sep="")) %>%
          left_join(coldf, by = "method") %>%
          left_join(shapedf, by = "method")
ddf[ddf$method=="rhoprop",]$association <- "propr"
ddf[ddf$method=="rhoshrink",]$association <- "corshrink"


## More helper functions for handling distance matrices
unpack <- function(x, dtypes) {
  if (inherits(x, 'list') && inherits(dtypes, 'list')) {
    return(do.call('rbind',
            lapply(1:length(x), function(i) unpack(x[[i]], dtypes[[i]]))
    ))
  }
  # x => distance matrix
  ex <- Matrix::Matrix(as.matrix(x), sparse=TRUE)
  es <- Matrix::summary(ex)
  es[,1] <- dtypes[es[,1]]
  es[,2] <- dtypes[es[,2]]
  ## filter identical data types
  es <- es[es[,1]==es[,2],]
  colnames(es) <- c("dtype", "dtype", "dist")
  as.data.frame(es[,2:3])
}
