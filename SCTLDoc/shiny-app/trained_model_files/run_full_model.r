wd <- getwd()

# Read in data
counts <- read.table('output_data/test_asv_counts.txt',header=TRUE,sep='\t')
metadata <- read.csv('input_data/SCTLD_meta_analysis_metadata.csv')

metadata <- metadata[!is.na(metadata$binary_disease),]
counts <- counts[,colnames(counts) %in% metadata$sample.id]
counts <- counts[apply(counts,1,sd) > 0,]

phy <- ape::keep.tip(ape::read.tree('output_data/reference_out_placement.tog.relabelled.tre'), rownames(counts))

logcounts <- log(as.matrix(counts))
logcounts[is.infinite(logcounts)] <- NA
size_factors <- apply(logcounts,2,function(x) sum(x,na.rm=TRUE) / length(x))
cl <- sapply(1:ncol(logcounts),function(x) logcounts[,x] - size_factors[[x]])

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
                 global_scale_prior = mean(apply(cl,1,sd,na.rm=T),na.rm=T),
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
                          'stepsize=0.001',
                          'engine=nuts',
                          'max_depth=8',
                          'num_warmup=1000',
                          'num_samples=1000',
                          'num_threads=4',
                          sep=' ')

setwd(file.path(wd,'models'))
print(sampling_command)
print(date())
system(sampling_command)
