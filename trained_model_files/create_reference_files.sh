## for now, use precomputed GreenGenes version of SEPP
cd output_data
## removed 'trap' command in run-sepp.sh so temp files are not removed
~/scripts/sepp-package/run-sepp.sh test_asv_seqs.fasta reference_out
