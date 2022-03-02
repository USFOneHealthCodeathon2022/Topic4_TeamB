data {
  int NF;
  int NS;
  array[NS] int disease;
  vector[NF] prior_intercepts;
  vector[NF] prior_disease_betas;
  matrix[NF,NF] prior_cov_intercepts;
  matrix[NF,NF] prior_cov_disease_betas;
  matrix[NF,NF] prior_cov_residuals;
  array[NS,NF] int counts;
}
transformed data {
  matrix[NF,NF] L_int = cholesky_decompose(prior_cov_intercepts);
  matrix[NF,NF] L_disease = cholesky_decompose(prior_cov_disease_betas);
  matrix[NF,NF] L_resid = cholesky_decompose(prior_cov_residuals);
}
parameters {
  vector[NF] intercepts;
  vector[NF] disease_betas;
  array[NS] vector[NF] latent_abundances;
}
model {
  intercepts ~ multi_normal_cholesky(prior_intercepts, L_int);
  disease_betas ~ multi_normal_cholesky(prior_disease_betas, L_disease);
  for(s in 1:NS) {
    latent_abundances[s] ~ multi_normal_cholesky(intercepts + disease[s] * disease_betas, L_resid);
    counts[s] ~ multinomial_logit(latent_abundances[s]);
  }
}
