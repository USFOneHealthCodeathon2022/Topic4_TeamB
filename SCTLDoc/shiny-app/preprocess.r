#Import packages
library(qiime2R)
library(tidyr)
library(tibble)
library(ggplot2)
library(phyloseq)

# Read qza files
#Import count table
ASV_table <- read_qza("input_data/tableV_BacArc_99_SCTLD.qza")
ASV_table <- ASV_table$data
ASV_table <- ASV_table[sample(1:nrow(ASV_table),1000),]
#Import taxon
taxa_16S <- read_qza("input_data/taxaVsearch_rep-seqs-dn-99_SCTLD.qza")
taxa_16S <- taxa_16S$data %>% as_tibble() %>% separate(Taxon, sep=";",
                                                       c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Convert taxonomy info to data frame with correct taxonomy labels
taxa_16S <- taxa_16S[-1,]

#Import metadata file
meta_data <- read.csv("input_data/SCTLD_meta_analysis_metadata.csv", header = T, row.names = 1,
                      na.strings = c("", "NA"))

#Make a phyloseq object
physeq_16S <- phyloseq(otu_table(ASV_table, taxa_are_rows= T),
                       tax_table(as.data.frame(taxa_16S) %>% column_to_rownames("Feature.ID") %>%
                                   as.matrix()), sample_data(meta_data))

ps <- subset_samples(physeq_16S, sample_type == "TissueSlurry" | sample_type == "Mucus" | sample_type == "TissueSlurry_Skleton" |
                     sample_type=="Seawater" | sample_type=="Sediment")
ps <- filter_taxa(ps, function(x) sum(x > 20) > (0.01*length(x)), TRUE)

save(ps, file='SCTLDoc/reference_data.RData')
