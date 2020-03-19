# Shrinkage improves estimation of microbial associations under different normalization methods

## Background

Consistent estimation of associations in microbial genomic survey count data is fundamental to microbiome research.
Technical limitations, including compositionality, low sample sizes, and technical variability, 
obstruct standard application of association measures and require data normalization prior to estimating associations.

In this project, we investigate the interplay between data normalization and microbial association estimation through
comprehensive analysis of statistical consistency. Leveraging the large sample size of the American Gut Project (AGP),
we assess the performance of the two linear association estimators, correlation and proportionality, under
different sample scenarios and data normalization schemes, including RNA-seq analysis work flows and log-ratio
transformations. 

![Project workflow](https://i.imgur.com/qYUXLy0.png)

We show that shrinkage estimation can universally improve the quality of association estimates for microbiome data and 
the consistency of downstream exploratory microbial data analysis tasks.

## Data and code availability 
This is the GitHub repository for the project. Please see below for all technical details.

The corresponding Synapse project is hosted at [syn21654780](https://www.synapse.org/#!Synapse:syn21654780). Data used for this study was accessed at ftp://ftp.microbio.me/AmericanGut/ag-2017-12-04/.

The latest complete American Gut Project dataset can be accessed on Qiita using study [ID 10317](https://qiita.ucsd.edu/study/description/10317).

## Repository structure
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


## Setup dependencies and data
Clone this repo and then set up conda environment with packages/dependencies
```sh
conda env create -f docker/environment.yml -n normcorr
conda activate normcorr
```

Or pull down and run via the pre-built docker container
```sh
docker pull docker.synapse.org/syn21654780/normcorr:latest
docker tag docker.synapse.org/syn21654780/normcorr:latest normcorr:latest
docker run -w $PWD -v $PWD:$PWD -ti normcorr:latest
```

### Sync data from Synapse
Data, code, and provenance for this project is hosted on [Synapse](synapse.org), project ID [syn21654780](https://www.synapse.org/#!Synapse:syn21654780). Get a synapse account to proceed.

The code is mirrored on GitHub, but the data is stored
on a public google drive.
Source files `data/amgut/` and intermediate RDS files `data/RDS/`
can be downloaded individually via the Synapse console or by
syncing via python scripts via the (experimental) [SynapseSync package](https://github.com/zdk123/SynapseSync) - included in the conda environment, or `pip install git+https://github.com/zdk123/SynapseSync.git` for the latest version.

Either set up a [synapseConfig file](https://python-docs.synapse.org/build/html/Credentials.html) or proceed with user name and password in a python session:

```python
from synapsesync import SynpaseProject, GDriveSession
from synapseutils.sync import syncFromSynapse

## With a ~/.synapseConfig file
syn = SynpaseProject("syn21654780")
## OR ##
## With user/pass
syn = SynpaseProject()
syn.login("<usr>", "<pass>")
syn.set_project("syn21654780")

syn.set_session(GDriveSession())

## sync source data and/or RDS files
syncFromSynapse(syn, "syn21654866", path="data/amgut/")
syncFromSynapse(syn, "syn21654865", path="data/RDS/")
```

## Reproduce analysis

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
```


## References

A prior version of the project has been posted on biorxiv:

[1] Michelle Badri, Zachary D. Kurtz, Christian L. Müller, Richard Bonneau, [Normalization methods for microbial abundance data strongly affect correlation estimates](https://www.biorxiv.org/content/10.1101/406264v1)

The current version of the project is available on biorxiv:

[2] Michelle Badri, Zachary D. Kurtz, Richard Bonneau, Christian L. Müller [Shrinkage improves estimation of microbial associations under different normalization methods]()


