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

## LFQ - reproducility 

Here we computing the reproducibility of LFQ scores within and between replicates

Usage:
```
/opt/R/4.2/bin/Rscript src/lfq_reproducibility.R ./data
```

Output:
```
figures/lfq_reproducibility.pdf
```


## SILAC + LFQ + SLiM-Domain Analysis

Here, we integrate SILAC and LFQ measurements in the context of SLiMs in screened peptides and cognate PFAM
domains in the interaction partners.

Usage:
```
/opt/R/4.2/bin/Rscript ./src/silac_lfq_slim_domain_analysis.R ./data/
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







