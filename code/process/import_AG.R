library(phyloseq)
library("doParallel")
registerDoParallel(cores=parallel::detectCores())


#setwd("/notrim/raw")

ag    <- import_biom("data/amgut/otu_table__BODY_HABITAT_UBERON_feces__.biom",parallel = TRUE)


map   <- read.delim("data/amgut/metadata__BODY_HABITAT_UBERON_feces__.txt", sep="\t", header=TRUE, row.names=1)
sample_data(ag) <- map

depths <- colSums(ag@otu_table@.Data)

ag.filt1 <- prune_samples(depths >= 10000, ag)

freq <- rowSums(sign(ag.filt1@otu_table@.Data))
ag.filt2 <- prune_taxa(freq > .3*nsamples(ag.filt1), ag.filt1)

depths   <- colSums(ag.filt2@otu_table@.Data)

ag.filt3 <- prune_samples(depths > quantile(depths, probs=seq(0, 1, .1))[2], ag.filt2)
depths   <- colSums(ag.filt3@otu_table@.Data)

X <- t(ag.filt3@otu_table@.Data)

ag.filt3
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 531 taxa and 9631 samples ]
# sample_data() Sample Data:       [ 9631 samples by 523 sample variables ]
# tax_table()   Taxonomy Table:    [ 531 taxa by 7 taxonomic ranks ]

saveRDS(ag.filt3, file="data/RDS/ag.filt3.RDS")

