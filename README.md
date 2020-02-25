## NormCorr
Compare normalization methods and correlation patterns for microbiome data


## Repo Structure:
```
.
├── code
│   ├── process
│   │   └── import_AG.R (import/filter American Gut scripts)
│   ├── correlation
│   │   ├── draw_samples.R (normalization and correlation methods)
│   │   ├── est_correlations.R (draw subsamples and estimate)
│   │   ├── shuffled_correlations.R (estimate null correlations)
│   │   └── plot_correlations.R (plot correlation distributions)
│   ├── compare
│   │   ├── compute_distances.R (parallelized correlation matrix distances)
│   │   ├── metrics.R (distance metrics in C++)
│   │   ├── data_helpers.R
│   │   ├── plot_distances.R (intra-method distances)
│   │   └── plot_embeddings.R (MDS plots)
│   ├── clustering
│   │   ├── spectral.R  (spectral clustering and tax compositions)
│   │   ├── hierarchical.R (hclust and circos plots)
│   │   └── plot_helpers.R
│   ├── networks
│   │   ├── make_nets.R (create relevance networks)
│   │   ├── plot_nets.R (igraph plotting)
│   │   ├── plot_lines.R (community analysis)
│   │   ├── consensus.R (upset plots, venn diagram and consensus network)
│   │   └── plot_helpers.R
│   ├── supp
│   │   ├── plot_msd.R (plot mean-sd plots)
│   │   ├── plot_msd.R (plot lambda values for cor shrinkage & rhoshrink)
│   │   └── plot_connectedcomponents.R (plot ranked eigenvalues)
│   └── helpers (helper functions and yaml files)
├── data
│   ├── amgut (American Gut project biom and mapping files)
│   └── RDS (intermediate data files)
├── plots (pdf figures)
└── docker (conda environment and container)
```


# Reproduce analysis

Set up conda environment with packages/dependencies
```sh
conda env create -f docker/environment.yml -n normcorr
conda activate normcorr
```
All scripts should be run from the base directory.

### Process American Gut data

Imports the biom file, filters and writes phyloseq RDS.
```sh
Rscript code/process/import_AG.R
```

### Estimate correlation/association measures
```
Rscript code/correlation/est_correlations.R
Rscript code/correlation/shuffled_correlations.R
Rscript code/correlation/plot_correlations.R
```

### Compare correlation matrices
```
Rscript code/compare/compute_distances.R
Rscript code/compare/plot_distances.R
Rscript code/compare/plot_embeddings.R
```

### Clustering
```
Rscript code/clustering/spectral.R
Rscript code/clustering/hierarchical.R
```

### Relevance networks
```
Rscript code/networks/make_nets.R
Rscript code/networks/plot_nets.R
Rscript code/networks/plot_lines.R
Rscript code/networks/consensus.R
```

### Supplemental figures
```
Rscript code/supp/plot_msd.R
Rscript code/supp/plot_lambda.R
Rscript code/supp/plot_connectedcomponents.R
```
