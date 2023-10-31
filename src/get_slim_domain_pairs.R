# given the peptide sequences and the array data which contains all detected peptide-protein interactions
# for all peptides and for all detected proteins in the array, find out which slim/domain pairs
# are compatible.
library(slimR)
library(data.table)
library(parallel)

datadir <- '/data/local/buyar/collaborations/trendelina/data'
peptides <- data.table::fread(file.path(datadir, 'from_trendelina', 'PeptideCandidates1.txt'))
peptides$type <- 'wt'
peptides[grepl('_', peptides$Name)]$type <- 'mut'
peptides[grepl('phos', peptides$Name),]$type <- 'phos'
# remove the leading lower-case 'p' in the phosphorylated form of the sequence
#peptides$Peptide <- gsub("p([STY])", paste0("_phos_", "\\1"), peptides$Peptide)
peptides$Peptide <- gsub("p([STY])", "\\1", peptides$Peptide)
# add leading and trailing "XXX" to the peptide sequence to avoid matching terminal motifs
peptides$Peptide <- paste0('XXX',peptides$Peptide,'XXX')
# assign experiment group ids to make it compatible with Expgroup ids in array data
peptides$Expgroup <- ceiling(peptides$ID/3)*3-2

# get elm classes
elms <- slimR::getElmClasses()

#get uniprot accessions and corresponding gene names from the fasta file of human uniprot sequences
message('Importing fasta sequences')
fasta <- Biostrings::readAAStringSet(
  slimR::getUniprotData(outDir = datadir, format = 'fasta', update = FALSE)
)

uniprot2geneName <- unique(data.table('uniprotAccession' = sub("^sp\\|(.+?)\\|.+$", "\\1", names(fasta)),
                                      'geneName' = sub("^.+ GN=(.+?) .+$", "\\1", names(fasta)),
                                      stringsAsFactors = FALSE))
# simplify fasta headers of fasta data
names(fasta) <- sub("^sp\\|(.+?)\\|.*", "\\1", names(fasta))

#import peptide-array results table from Katrina
arrayDataFile <- file.path(datadir, 'from_trendelina', 'SILACLFQ_extended_Bora_loose_p001FD_1.txt')
arrayData <- fread(arrayDataFile)
# simplify uniprot accession/gene name fields => remove unreviewed ids and keep only one per row
cl <- parallel::makeCluster(10)
parallel::clusterExport(cl, varlist = c('arrayData', 'uniprot2geneName'))
s <- pbapply::pbsapply(cl = cl, strsplit(arrayData$Majority.protein.IDs, ";"), function(x) {
  # pick first match
  intersect(x, uniprot2geneName$uniprotAccession)[1]
})
arrayData$uniprotAccession <- s
parallel::stopCluster(cl)

# modify arrayData to have one line per protein & peptide in the screen
arrayData$SILACGroup <- gsub("[0-9]+$", "", arrayData$SILACGroup)
colnames(arrayData) <- gsub("_", ".", gsub("LFQsignificant", "lfq_", colnames(arrayData)))

arrayData <- dcast.data.table(arrayData[!is.na(uniprotAccession)],
                 Expgroup + uniprotAccession + lfq.Wt.loose + lfq.Phos.loose + lfq.Mut.loose ~ SILACGroup,
                 value.var = 'Median.SILAC.ratio')
# keep interactions that are not NA for any silac values
to_keep <- apply(arrayData[,c('phos_mut', 'wt_mut', 'wt_phos')], 1, function(x) sum(!is.na(x)) == 3)
arrayData <- arrayData[to_keep]

message("importing PFAM domains")
uniprot2pfam <- as.data.table(slimR::getPFAM(organism = 9606, pfam_version = "Pfam34.0"))
uniprot2pfam <- uniprot2pfam[uniprot2pfam$seqnames %in% unique(arrayData$uniprotAccession),]
uniprot2pfam <- unique(subset(uniprot2pfam, select = c('seqnames', 'pfam_acc', 'pfam_name', 'clan')))
colnames(uniprot2pfam) <- c('uniprotAccession', 'PFAM', 'NAME', 'CLAN')
pfamClans <- slimR::getPFAMClans()
message("importing ELM -> PFAM associations")
# download a table of ELM classes and their cognate PFAM domains
elm2pfam <- data.table(slimR::getELMdomainInteractions())

# simplify uniprot2pfam table:
# remove proteins that don't exist in the screen
uniprot2pfam <- uniprot2pfam[uniprotAccession %in% arrayData$uniprotAccession]


# for each peptide, and for each protein in the screen, find which elm-interacting domains exist
slimDomainInteractions <- do.call(rbind, pbapply::pblapply(unique(peptides$Expgroup), function(e) { #for each experiment group
  # get list of proteins detected for the experiment group
  prots <- unique(arrayData[Expgroup == e]$uniprotAccession)
  dt <- do.call(rbind, lapply(c(unique(peptides$type)), function(g) { # for each genotype
    # slims in peptide
    seq <- peptides[Expgroup == e][type == g]$Peptide
    slims <-sort(unique(slimR::searchSLiMs(seq, elms)[['SLiM']]))
    # find elm-interacting domains among the detected interactors
    dt <- droplevels(merge(uniprot2pfam[uniprotAccession %in% prots],
          elm2pfam[ELM_identifier %in% slims], by.x = 'NAME', by.y = 'Interaction_Domain_Name'))
    # get a list of slim-domain pairs for each protein
     do.call(rbind, lapply(unique(dt$uniprotAccession), function(x) {
      data.table('uniprotAccession' = x,
                 'slim_domain_pairs' =  paste(paste(dt[uniprotAccession == x]$NAME,
                                                    dt[uniprotAccession == x]$ELM_identifier, sep = "->"),
                                              collapse = '; '),
                 'genotype' = paste0(g, '.slim_domain_pairs'),
                 'Expgroup' = e)
    }))
  }))
}))

slimDomainInteractions <- dcast.data.table(slimDomainInteractions,
                                           Expgroup + uniprotAccession ~ genotype,
                                           value.var = 'slim_domain_pairs')

dt <- merge(arrayData, slimDomainInteractions, by = c('Expgroup', 'uniprotAccession'))

# find slim-domain pairs gained/lost in mutant peptide compared to wt/phos
l1 <- strsplit(dt$wt.slim_domain_pairs, '; ')
l2 <- strsplit(dt$mut.slim_domain_pairs, '; ')
dt <- cbind(dt,do.call(rbind, lapply(1:length(l1), function(i) {
  gained <- setdiff(na.omit(l2[[i]]), na.omit(l1[[i]]))
  lost <- setdiff(na.omit(l1[[i]]), na.omit(l2[[i]]))
  gained <- ifelse(length(gained) > 0, gained, NA)
  lost <- ifelse(length(lost) > 0, lost, NA)
  data.table('gained_in_mut'= gained, 'lost_in_mut' = lost)
})))

library(ggpubr)
library(ggplot2)
ggplot2::theme_set(ggpubr::theme_pubclean())


p1 <- ggpubr::ggboxplot(dt, y = 'wt_phos', add = 'jitter', color = 'lfq.Phos.loose') +
  labs(title = 'Interactions which can be explained \n by SLiM-Domain pairs (WT vs Phos)',
       y = 'SILAC: log2(WT/Phos)') +
  geom_hline(yintercept = 0, size = 4, alpha = 0.2, color = 'black') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_color_brewer(type = 'qual', palette = 3, direction = -1) +
  ylim(-6, 6)

p2 <- ggpubr::ggboxplot(arrayData, y = 'wt_phos', add = 'jitter', color = 'lfq.Phos.loose') +
  labs(title = 'All interactions (WT vs Phos)',
       y = 'SILAC: log2(WT/Phos)') +
  geom_hline(yintercept = 0, size = 4, alpha = 0.2, color = 'black') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_color_brewer(type = 'qual', palette = 3, direction = -1) +
  ylim(-6, 6)


p <- cowplot::plot_grid(p1, p2, align = 'h')

ggsave(filename = 'slimDomain.wt_vs_phos.pdf', plot = p,
       units = 'in', width = 10, height = 8)

p3 <- ggpubr::ggboxplot(dt, y = 'phos_mut', add = 'jitter', color = 'lfq.Phos.loose') +
  labs(title = 'Interactions which can be explained \n by SLiM-Domain pairs (Phos vs Mut)',
       y = 'SILAC: log2(Phos/Mut)') +
  geom_hline(yintercept = 0, size = 4, alpha = 0.2, color = 'black') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_color_brewer(type = 'qual', palette = 3, direction = -1) + 
  ylim(-6, 6)

p4 <- ggpubr::ggboxplot(arrayData, y = 'phos_mut', add = 'jitter', color = 'lfq.Phos.loose') +
  labs(title = 'All interactions (Phos vs Mut)',
       y = 'SILAC: log2(Phos/Mut)') +
  geom_hline(yintercept = 0, size = 4, alpha = 0.2, color = 'black') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_color_brewer(type = 'qual', palette = 3, direction = -1) +
  ylim(-6, 6)


p <- cowplot::plot_grid(p3, p4)
ggsave(filename = 'slimDomain.phos_vs_mut_1.pdf', plot = p,
       units = 'in', width = 10, height = 8)


dt$mutation_breaks_slimDomain_pair <- !is.na(dt$lost_in_mut)
p5 <- ggpubr::ggboxplot(dt, y = 'phos_mut', add = 'jitter', facet.by = 'lfq.Phos.loose',
                  color = 'mutation_breaks_slimDomain_pair') +
  labs(title = 'Interactions which can be explained \n by SLiM-Domain pairs (Phos vs Mut)',
       y = 'SILAC: log2(Phos/Mut)') +
  geom_hline(yintercept = 0, size = 4, alpha = 0.2, color = 'black') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_color_brewer(type = 'qual', palette = 1)

ggsave(filename = 'slimDomain.phos_vs_mut_2.pdf', plot = p5,
       units = 'in', width = 10, height = 8)


# # import localisation data
# loc.dt <- data.table::fread(file.path(datadir, 'deeploc', 'deeploc.tsv'))
# arrayData$localisation <- loc.dt[match(arrayData$uniprotAccession, uniprotAccession)]$localisation.simple
# arrayData[is.na(localisation)]$localisation <- 'unknown'
#
# dt <- arrayData[,length(unique(uniprotAccession)),by = c('Expgroup', 'localisation', 'lfq.Phos.loose')]
#
# dtc <- dcast(dt, Expgroup + lfq.Phos.loose ~ localisation, value.var = 'V1', fill = 0)
#
# # ggplot(dt, aes(y = wt_phos)) + geom_boxplot(aes(color = lfq.Phos.loose))
# # ggplot(arrayData, aes(y = wt_phos)) + geom_boxplot(aes(color = lfq.Phos.loose))
# # ggplot(arrayData, aes(y = wt_phos)) + geom_boxplot(aes(color = lfq.Wt.loose))
# #
# # ggplot(dt, aes(y = phos_mut)) + geom_boxplot(aes(color = lfq.Phos.loose))
# # ggplot(dt[!is.na(gained_in_mut)], aes(y = phos_mut)) +
# #   geom_boxplot(aes(color = lfq.Mut.loose))
# # ggplot(dt[!is.na(lost_in_mut)], aes(y = phos_mut)) +
# #   geom_boxplot(aes(color = lfq.Phos.loose))
#
#
#
#
# # explore protein embeddings for interaction partners detected in the array
# # e.g. can we predict from sequence embeddigns if the proteins will bind preferentially to
# # the wt/mut/phos forms of the peptide?
# source('/data/local/buyar/arcas/uyar_et_al_multiomics_deeplearning/src/common_functions.R')
# Embeddings <- data.table::fread('/data/local/buyar/arcas/protein_embeddings/prot_trans/embeddings/uniprot/embeddings.protein_level.csv')
# arrayData$lfq.cat <- paste0("mut", arrayData$lfq.Mut.loose,
#                             "phos", arrayData$lfq.Phos.loose,
#                             "wt", arrayData$lfq.Wt.loose)
# arrayData[,length(unique(uniprotAccession)),by = c('lfq.Phos.loose', 'Expgroup')][order(V1, decreasing = T)][lfq.Phos.loose == '+'][1:10]
# e <- 13#22#103 #13
# selected <- unique(arrayData[Expgroup == e][lfq.cat != 'mut-phos-wt-']$uniprotAccession)
# E <- as.matrix(subset(Embeddings, select = selected))
# dt <- arrayData[Expgroup == e][lfq.cat != 'mut-phos-wt-'][match(colnames(E), uniprotAccession)]
# df <- plot_tsne(t(E))
# df <- cbind(df, dt[,c('phos_mut', 'wt_mut', 'wt_phos')])
#
# ggplot(df, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = scale(wt_phos) > 0))
# ggplot(df, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = scale(wt_mut) > 0))
# ggplot(df, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = scale(phos_mut) > 0))
#
# rocs <- predict_cluster_labels_glmnet(t(E), as.factor(f), returnROC = T, repeats = 5, sampling = 'down')
# plot_roc_list(rocs)
#
# # other interesting questions:
# # 1. Is there a difference in disorder content of proteins that preferentially bind different forms of peptides
#
#
#
#
#
#
#
#
#
#
#
#
