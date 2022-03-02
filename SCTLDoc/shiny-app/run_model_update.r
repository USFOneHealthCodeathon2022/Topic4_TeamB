wd <- getwd()

args <- commandArgs(TRUE)

# Read in data
counts <- read.table(args[[1]],header=TRUE,sep='\t')
metadata <- read.csv(args[[2]])
full_phy <- ape::read.tree(args[[3]])

ref_tips <- ''
new_tips <- rownames(counts)
phy_intersect <- ape::keep.tip(full_phy,c(ref_tips,new_tips))

prior_corr_intersect <- ape::vcv(phy_intersect,corr=TRUE)

reference_prior_corr <- read.table(file.path('trained_model_files','input_correlation.txt'), header=TRUE, sep='\t')
reference_prior_corr_L <- t(chol(reference_prior_corr))
stan.fit <- cmdstanr::read_cmdstan_csv(Sys.glob(path.expand(file.path('trained_model_files', 'full_model_samples_*.csv'))), format = 'draws_array')

idx_intercepts <- grep('intercepts',dimnames(stan.fit$post_warmup_draws)[[3]])
idx_disease_betas <- grep('disease_betas',dimnames(stan.fit$post_warmup_draws)[[3]])
idx_residuals <- grep('residuals',dimnames(stan.fit$post_warmup_draws)[[3]])

prior_intercepts <- apply(stan.fit$post_warmup_draws[,,idx_intercepts], 3, mean)
prior_disease_betas <- apply(stan.fit$post_warmup_draws[,,idx_disease_betas], 3, mean)
prior_residuals <- apply(stan.fit$post_warmup_draws[,,idx_residuals], 3, mean)

old_dims <- dim(stan.fit$post_warmup_draws)
old_NF <- length(idx_intercepts)
new_NF <- length(new_tips)

reference_vs_new_prior_corr <- prior_corr_intersect[ref_tips,new_tips]
new_prior_corr <- prior_corr_intersect[new_tips,new_tips]
new_intercept_samples <- array(NA,c(old_dims[[1]],old_dims[[2]],new_NF))
new_disease_betas_samples <- array(NA,c(old_dims[[1]],old_dims[[2]],new_NF))
new_disease_betas_samples <- array(NA,c(old_dims[[1]],old_dims[[2]],new_NF))
## https://betanalpha.github.io/assets/case_studies/gp_part1/part1.html
for(chain in 1:old_dims[[2]]) {
  for(iter in 1:old_dims[[1]]) {
    intercept_scale <- stan.fit$post_warmup_draws[iter,chain,'intercept_scale']^2
    reference_prior_cov_L <- intercept_scale * reference_prior_corr_L
    K_div_intercepts_old <- forwardsolve(reference_prior_cov_L, stan.fit$post_warmup_draws[iter,chain,idx_intercepts])
    K_div_intercepts_old <- backsolve(reference_prior_cov_L, K_div_intercepts_old, transpose=TRUE)
    reference_vs_new_prior_cov_intercept <- intercept_var * reference_vs_new_prior_corr
    f2_mu <- K_div_intercepts_old %*% reference_vs_new_prior_cov_intercept
    v_pred <- forwardsolve(reference_prior_cov_L, reference_vs_new_prior_cov_intercept)
    cov_f2 <- intercept_var * new_prior_corr - tcrossprod(v_pred)
    new_intercept_samples[iter,chain,] <- rnorm(new_NF) %*% chol(cov_f2)
  }
}


prior_cov_intercepts <- matrix(0,nrow=old_NF,ncol=old_NF)
prior_cov_disease_betas <- matrix(0,nrow=old_NF,ncol=old_NF)
prior_cov_residuals <- matrix(0,nrow=old_NF,ncol=old_NF)
for(chain in 1:old_dims[[2]]) {
  for(iter in 1:old_dims[[1]]) {
    diff <- stan.fit$post_warmup_draws[iter,chain,idx_intercepts] - prior_intercepts
    prior_cov_intercepts <- prior_cov_intercepts + tcrossprod(diff)
    diff <- stan.fit$post_warmup_draws[iter,chain,idx_disease_betas] - prior_disease_betas
    prior_cov_disease_betas <- prior_cov_disease_betas + tcrossprod(diff)
    diff <- stan.fit$post_warmup_draws[iter,chain,idx_residuals] - prior_residuals
    prior_cov_residuals <- prior_cov_residuals + tcrossprod(diff)
  }
}
prior_cov_intercepts <- prior_cov_intercepts / (prod(old_dims[1:2])-1)
prior_cov_disease_betas <- prior_cov_disease_betas / (prod(old_dims[1:2])-1)
prior_cov_residuals <- prior_cov_residuals / (prod(old_dims[1:2])-1)


relabund <- apply(counts,2,function(x) x/sum(x))
clr <- apply(relabund,2,function(x) log(x) - mean(log(x[x>0])))
clr[is.infinite(clr)] <- NA

disease <- metadata$binary_disease
names(disease) <- metadata$sample.id
disease <- disease[!is.na(disease)]
disease[disease %in% c('n','N')] <- 0
disease[disease %in% c('y','Y')] <- 1
disease <- as.numeric(disease)

counts <- counts[,colnames(counts) %in% names(disease)]
disease <- disease[colnames(counts)]

prior_corr <- ape::vcv(phy,corr=TRUE)

stan_dat <- list(NF=nrow(counts),
                 NS=ncol(counts),
                 disease=disease,
                 global_scale_prior = mean(apply(clr,1,sd,na.rm=T),na.rm=T),
                 prior_corr = prior_corr,
                 counts = t(counts))

cmdstanr::write_stan_json(stan_dat, file.path('trained_model_files', 'full_model_input.json'))

setwd(cmdstanr::cmdstan_path())
system(paste0('make ', file.path(wd,'models','full_model')))

sampling_command <- paste('./full_model',
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
