[![DOI](https://zenodo.org/badge/712373714.svg)](https://zenodo.org/doi/10.5281/zenodo.10786077)

# Description

Data and scripts related to short linear motif analyses for the collaboration project titled 
**"Pathogenic mutations of human phosphorylation sites affect protein-protein interactions"** with Trendelina Rrustemi from Matthias Selbach's Lab.

# Publications

See the manuscript on Biorxiv [here](https://www.biorxiv.org/content/10.1101/2023.08.01.551433v1.full).

# Installation

Clone the repo: 
```
git clone git@github.com:BIMSBbioinfo/collab_rrustemi_selbach_prisma.git
```

# Managing Dependencies

To run the scripts within this repo, you need an `R (version >= 4.2)` installation. 

## Using `renv` package 

You can create an environment using the environment snapshot file `renv.lock` in the existing repo folder. 
However, this requires the `renv` packaged to be installed. 

```
# create an R session
R
# install the renv package
install.packages('renv')
```

Once, `renv` is installed, you can use the `renv.lock` file in the repo to restore the snapshot of
the environment with all the necessary packages. 

```
Rscript -e "library(renv); renv::init(); renv::restore()"
```

To deactivate the session
```
Rscript -e "renv::deactivate()"
```


## Manual Installation

Alternatively, the dependencies can be installed using `BiocManager` and `devtools` packages.  

### CRAN and Bioconductor Packages

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('ggplot2', 'data.table', 'ggpubr', 'ComplexHeatmap', 'cowplot', 'parallel', 'GenomicRanges', 'Biostrings', 'rmarkdown', 'knitr', 'pbapply'))

```

### Devtools 

```
install.packages('devtools') 
devtools::install_github('BIMSBbioinfo/slimR')
```

# Analyses

## LFQ - reproducility 

Here we computing the reproducibility of LFQ scores within and between replicates

Usage:
```
Rscript src/lfq_reproducibility.R ./data `pwd` 
```

Output:
```
figures/lfq_reproducibility.pdf
```

## LFQ + SLiM Domain Analysis 

The goal of this analysis is to inspect the LFQ scores in the context of SLiM-Domain
interactions. 

This analysis is done within an rmarkdown file which includes code, text, and figures.

Usage:
```
Rscript -e "rmarkdown::render('src/LFQ_slim_domain_analysis.Rmd', output_dir = './figures')"
```

Output:
```
- figures/LFQ_slim_domain_analysis.html
- figures/lfq_zscore_phos_vs_wt_vs_mut.doc_lig.pdf  
- figures/lfq_zscore_phos_vs_wt_vs_mut.doc_lig_deg.pdf  
- figures/lfq_zscore_vs_slimdomain_interactions.pdf
- tables/LFQinteractions.scaled_by_peptide.tsv  
- tables/LFQscaled_slim_domain_interactions.tsv

```

## Phospho-dependent Domain Enrichment Analysis

This analysis aims to observe if we can detect an enrichment of phospho-dependent binding domains in the array data
among proteins that are preferentially binding to phosphorylated forms of the peptides. 

Usage: 
```
Rscript -e "rmarkdown::render('src/phospho_domain_discovery.Rmd', output_dir = './figures')"
```

Output:
```
See figures/phospho_domain_analysis.pdf and figures/phospho_domain_discovery.html
```




## SILAC + LFQ + SLiM-Domain Analysis

Here, we integrate SILAC and LFQ measurements in the context of SLiMs in screened peptides and cognate PFAM
domains in the interaction partners.

Usage:
```
Rscript ./src/silac_lfq_slim_domain_analysis.R ./data/ `pwd`
```  

We looked for slims in the peptides and PFAM domains found in the interaction partners that can 
bind those slims from the array data. If we subset the array data by if the interaction can be explained 
by such slim-domain pairs, then we can see significant differences silac ratio distributions 
when comparing wt vs phos, and phos vs mut. 

```
See figures/slimDomain.wt_vs_phos.pdf and figures/slimDomain.phos_vs_mut_1.pdf
```

Another interesting observation is that if we break down these slim-domain pairs into further groups 
such as if the mutant peptide loses a slim-domain pair (as the mutation breaks the motif pattern), 
then we see an even increased difference in silac ratios for phos vs mut peptides. 
```
See figures/slimDomain.phos_vs_mut_2.pdf
```










