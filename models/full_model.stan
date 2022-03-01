data {
  int NF;
  int NS;
  array[NS] int disease;
  real global_scale_prior;
  matrix[NF,NF] prior_corr;
  array[NS,NF] int counts;
}
transformed data {
  matrix[NF,NF] L = cholesky_decompose(prior_corr);
}
parameters {
  real<lower=0> global_scale;
  real<lower=0> intercept_scale;
  real<lower=0> disease_scale;
  vector[NF] intercepts;
  vector[NF] disease_betas;
  array[NS] vector[NF] latent_abundances;
}
model {
  matrix[NF,NF] L_resid = global_scale * L;
  global_scale ~ normal(0, global_scale_prior);
  intercept_scale ~ normal(0, 10*global_scale);
  disease_scale ~ normal(0, global_scale);
  intercepts ~ multi_normal_cholesky(rep_vector(0,NF), intercept_scale * L);
  disease_betas ~ multi_normal_cholesky(rep_vector(0,NF), disease_scale * L);
  for(s in 1:NS) {
    latent_abundances[s] ~ multi_normal_cholesky(intercepts + disease[s] * disease_betas, L_resid);
    counts[s] ~ multinomial_logit(latent_abundances[s]);
  }
}
generated quantities {
  matrix[NF,NS] residuals;
  for(s in 1:NS) {
    residuals[,s] = latent_abundances[s] - intercepts + disease[s] * disease_betas;
  }
}
