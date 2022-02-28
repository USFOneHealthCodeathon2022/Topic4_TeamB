args <- commandArgs(TRUE)

# Read in the 16S qza files
counts <- read.table(args[[1]],header=TRUE,sep='\t')
metadata <- read.csv('output_data/SCTLD_meta_analysis_metadata.csv')
phy <- phytools::read_newick('output_data/phy.newick')


