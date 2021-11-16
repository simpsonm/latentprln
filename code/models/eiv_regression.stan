data {
  int<lower = 1> nmetro;
  int<lower = 1> nx;
  vector[nmetro] index;
  vector<lower = 0>[nmetro] index_se;
  matrix[nmetro, nx] x;
  matrix<lower = 0>[nmetro, nx] x_se;
  real alpha_cs_prior_mean;
  real<lower = 0> alpha_cs_prior_sd;
  vector[nx] beta_cs_prior_mean;
  vector<lower = 0>[nx] beta_cs_prior_sd;
  vector[nx] mu_cs_prior_mean;
  vector<lower = 0>[nx] mu_cs_prior_sd;
  real tau_cs_prior_mean;
  real<lower = 0> tau_cs_prior_sd;
  vector[nx] sigma_cs_prior_mean;
  vector<lower = 0>[nx] sigma_cs_prior_sd;
}
transformed data {
  vector[nmetro] index_cs;
  vector<lower = 0>[nmetro] index_cs_se;
  matrix[nmetro, nx] x_cs;
  matrix<lower = 0>[nmetro, nx] x_cs_se;
  real mean_index;
  real sd_index;
  vector[nx] mean_x;
  vector[nx] sd_x;
  matrix<lower = 0>[nmetro, nx] x_var;
  vector<lower = 0>[nmetro] index_var;

  mean_index = mean(index);
  sd_index = sd(index);
  index_cs = (index - mean_index) / sd_index;
  index_cs_se = index_se / sd_index;

  for(n in 1:nx){
    mean_x[n] = mean(x[,n]);
    sd_x[n] = sd(x[,n]);
    x_cs[,n] = (x[,n] - mean_x[n]) / sd_x[n];
    x_cs_se[,n] = x_se[,n] / sd_x[n];
  }

  x_var = x_cs_se .* x_cs_se;
  index_var = index_cs_se .* index_cs_se;
}
parameters {
  real alpha_cs;
  vector[nx] beta_cs;
  vector[nx] mu_cs;
  vector<lower = 0>[nx] sigma_cs;
  real<lower = 0> tau_cs;  
}
transformed parameters {
  vector<lower = 0>[nx] sigma2;
  real<lower = 0> tau2;
  tau2 = tau_cs * tau_cs;
  sigma2 = sigma_cs .* sigma_cs;
}
model {
  for(i in 1:nx){
    x_cs[,i] ~ normal(mu_cs[i], sqrt(sigma2[i] + x_var[,i]));
  }
  for(i in 1:nmetro){
    real mubar;
    real sig2bar;
    vector[nx] lambda;
    lambda = sigma2 ./ (sigma2 + to_vector(x_var[i]) );
    mubar = alpha_cs + dot_product(mu_cs + lambda .* (to_vector(x[i]) - mu_cs), beta_cs);
    sig2bar = index_var[i] + tau2 +
      dot_product(sigma2 .* (1 - lambda), beta_cs .* beta_cs);
    index_cs ~ normal(mubar, sqrt(sig2bar));
  }
}
generated quantities {
  real alpha;
  vector[nx] beta;
  vector[nx] mu;
  vector<lower = 0>[nx] sigma;
  real<lower = 0> tau;

  for(n in 1:nx){
    beta[n] = (sd_index / sd_x[n]) * beta_cs[n];
    mu[n] = mu_cs[n] * sd_x[n] + mean_x[n];
    sigma[n] = sd_x[n] * sigma_cs[n];
  }
  tau = sd_index * tau_cs;
  alpha = sd_index * alpha_cs + sd_index - mean_x'*beta;
}
