data {
  int NF;
  int NS;
  row_vector[NS] disease;
  real global_scale_prior;
  real ou_prior;
  matrix[NF,NF] phy_dists;
  int N_nodes;
  array[N_nodes,2] int idx_descendants;
  array[N_nodes] int idx_self;
  array[NS,NF] int counts;
}
parameters {
  real<lower=0> gs_raw;
  real<lower=0> is_raw;
  real<lower=0> ds_raw;
  vector[NF] i_raw;
  vector[N_nodes] i_bal_raw;
  vector[N_nodes] bd_bal_raw;
  vector[N_nodes] res_bal_raw;
  vector[NF] db_raw;
  matrix[NF,NS] res_raw;
  real<lower=0> ou;
}
transformed parameters {
  real<lower=0> global_scale = global_scale_prior * gs_raw;
  real<lower=0> intercept_scale = global_scale * is_raw;
  real<lower=0> disease_scale = global_scale * ds_raw;
  matrix[NF,NF] L = cholesky_decompose(add_diag(exp(-phy_dists / ou), 1e-9));
  vector[NF] intercepts = intercept_scale * (L * i_raw);
  vector[NF] disease_betas = disease_scale * (L * db_raw);
  matrix[NF,NS] residuals = global_scale * (L * res_raw);
  matrix[NF,NS] balances = rep_matrix(i_bal_raw,NS) + bd_bal_raw * disease + res_bal_raw;
  for(s in 1:NS) {
    for(n in 1:N_nodes) {
      // start at tips and work toward root
      // not quite right because higher nodes have deterministic effects
      balances[idx_descendants[n,1]:idx_descendants[n,2],s] = balances[idx_descendants[n,1]:idx_descendants[n,2],s] - log_sum_exp(balances[idx_descendants[n,1]:idx_descendants[n,2],s]);
    }
  }
}
model {
  matrix[NF,NS] latent_abundances = rep_matrix(intercepts,NS) + balances + disease_betas * disease + residuals;
  gs_raw ~ std_normal();
  is_raw ~ std_normal();
  ds_raw ~ std_normal();
  i_raw ~ std_normal();
  db_raw ~ std_normal();
  to_vector(res_raw) ~ std_normal();
  ou ~ exponential(1);
  for(s in 1:NS) {
    counts[s] ~ multinomial_logit(latent_abundances[,s]);
  }
}
