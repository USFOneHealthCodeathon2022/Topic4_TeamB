# Read in the 16S qza files
ASVtable_16S <- qiime2R::read_qza("input_data/tableV_BacArc_99_SCTLD.qza")$data
seqs <- qiime2R::read_qza('input_data/rep-seqs-dn-99.qza')$data
rownames(ASVtable_16S) <- as.character(seqs[match(rownames(ASVtable_16S), names(seqs))])
names(seqs) <- as.character(seqs)

rel <- apply(ASVtable_16S,2,function(x) x/sum(x))
test_set <- ASVtable_16S[order(rowSums(rel),decreasing=TRUE)[1:1000],]
seqs <- seqs[rownames(test_set)]
write.table(test_set,file='output_data/test_asv_counts.txt', sep='\t', quote=FALSE)
Biostrings::writeXStringSet(seqs,'output_data/test_asv_seqs.fasta')
