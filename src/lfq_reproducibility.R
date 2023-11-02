# assess LFQ reproducibility within between replicates

library(ggplot2)
library(data.table)
library(ggpubr)

args <- commandArgs(trailingOnly = T)

datadir <- args[1]
outdir <- args[2]

if(is.na(datadir)) {stop("Provide a path to data folder")}
if(is.na(outdir)) {stop("Provide a path to output folder")}
if(!dir.exists(datadir)) {stop(datadir, "doesn't exist")}
if(!dir.exists(outdir)) {stop(outdir, "doesn't exist")}

data.table::setDTthreads(12)

# import LFQ data 
lfqData <- data.table::fread(file.path(datadir, 'LFQinteractions.txt'))

# compute pairwise correlations 
M <- as.matrix(lfqData[,-c(1:3)]) # remove non-numeric columns 
Mc <- cor(M, use = 'complete.obs')

# reshape and process 
dt <- data.table::data.table(data.table::melt(as.data.table(Mc, keep.rownames = T)))
colnames(dt)[1:2] <- c('Var1', 'Var2')
dt <- dt[Var1 != Var2]

dt$comp <- apply(dt, 1, function(x) paste(sort(x[1:2]), collapse = '_'))
# remove double counts
dt2 <- dt[,.SD[which.min(value)],by = comp]

# check if compared pulldowns are replicates
dt2$replicates <- gsub("\\.[12]$", "", dt2$Var1) == dt2$Var2
dt2$replicates <- ifelse(dt2$replicates == TRUE, "Replicate Pairs", "Non-replicate Pairs")

# make plot and save
p <- ggviolin(dt2, x = 'replicates', y = 'value', fill = 'replicates') + 
  stat_compare_means(label.x = 1.4, label.y = 0.8) + 
  labs(x = "Peptide pairs", y = "Pearson's correlation", fill = '') +
  scale_fill_brewer(type = 'qual', palette = 3) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))

figuredir <- file.path(outdir, "figures")
if (!dir.exists(figuredir)) { dir.create(figuredir) }

ggsave(plot = p, file.path(figuredir, "lfq_reproducibility.pdf"), 
       width = 9.42, height = 8.4, units = 'in')
