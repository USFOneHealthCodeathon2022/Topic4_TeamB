# Read in the 16S qza files
ASVtable_16S <- qiime2R::read_qza("input_data/tableV_BacArc_99_SCTLD.qza")$data
seqs <- qiime2R::read_qza('input_data/rep-seqs-dn-99.qza')$data
rownames(ASVtable_16S) <- as.character(seqs[match(rownames(ASVtable_16S), names(seqs))])

## add sequences to reference tree
system(paste0('python ~/scripts/sepp/run_sepp.py -t ', ref_tree, ' -r ', raxml_file, ' -a ', alignment_reference, ' -f ', new_seqs))



library(tidyr)
ASVtaxa_16S <- qiime2R::read_qza("input_data/taxaVsearch_rep-seqs-dn-99_SCTLD.qza")$data
taxtable_16S <- ASVtaxa_16S %>% as_tibble() %>% separate(Taxon, sep=";",c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) 
taxtable_16S <- as.data.frame(taxtable_16S)


# Convert taxonomy info to data frame with correct taxonomy labels
taxtable_16S <- taxtable_16S[-1,]
#asv_tab$asv_id <- rownames(asv_tab) # add a new column for ids

tree_file <-read_tree("input_data/exported-tree//tree.nwk")
