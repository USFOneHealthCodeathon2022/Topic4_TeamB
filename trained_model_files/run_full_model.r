wd <- getwd()

# Read in data
counts <- read.table('output_data/test_asv_counts.txt',header=TRUE,sep='\t')
metadata <- read.csv('input_data/SCTLD_meta_analysis_metadata.csv')
phy <- ape::keep.tip(ape::read.tree('output_data/reference_out_placement.tog.relabelled.tre'), rownames(counts))

metadata <- metadata[!is.na(metadata$binary_disease),]
counts <- counts[,colnames(counts) %in% metadata$sample.id]

relabund <- apply(counts,2,function(x) x/sum(x))
clr <- apply(relabund,2,function(x) log(x) - mean(log(x[x>0])))
clr[is.infinite(clr)] <- NA

disease <- metadata$binary_disease
names(disease) <- metadata$sample.id
disease[disease %in% c('n','N')] <- 0
disease[disease %in% c('y','Y')] <- 1
disease <- disease[colnames(counts)]
disease <- as.numeric(disease)


phy_dists <- ape::cophenetic.phylo(phy)
phy_dists <- phy_dists / max(phy_dists)

stan_dat <- list(NF=nrow(counts),
                 NS=ncol(counts),
                 disease=disease,
                 global_scale_prior = mean(apply(clr,1,sd,na.rm=T),na.rm=T),
                 ou_prior = mean(phy_dists),
                 phy_dists = phy_dists,
                 counts = t(counts))

cmdstanr::write_stan_json(stan_dat, file.path('trained_model_files', 'full_model_input.json'))

setwd(cmdstanr::cmdstan_path())
system(paste0('make ', file.path(wd,'models','full_model_nc')))

sampling_command <- paste('./full_model_nc',
                          paste0('data file=',path.expand(file.path('../trained_model_files', 'full_model_input.json'))),
                          'output',
                          paste0('file=',path.expand(file.path('../trained_model_files', 'full_model_samples.csv'))),
                          paste0('refresh=', 1),
                          'method=sample',
                          'num_chains=4',
                          'algorithm=hmc',
                          'engine=nuts',
                          'num_warmup=1000',
                          'num_samples=1000',
                          'num_threads=4',
                          sep=' ')

setwd(file.path(wd,'models'))
print(sampling_command)
print(date())
system(sampling_command)
