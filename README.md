# collab_rrustemi_selbach_prisma
Data and scripts related to short linear motif analyses for the collaboration with Trendelina Rrustemi from Matthias Selbach's Lab.

# Dependencies

## CRAN and Bioconductor Packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('ggplot2', 'data.table', 'ggpubr')

## Devtools 

```
install.packages('devtools') 
devtools::install_github('BIMSBbioinfo/slimR')
```

# Analyses

## Computing the reproducibility of LFQ scores within and between replicates

Usage:

```
/opt/R/4.2/bin/Rscript src/lfq_reproducibility.R ./data
```






