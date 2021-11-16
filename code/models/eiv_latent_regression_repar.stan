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
}
parameters {
  real alpha_cs;
  vector[nx] beta_cs;
  vector[nx] mu_cs;
  vector<lower = 0>[nx] sigma_cs;
  real<lower = 0> tau_cs;
  matrix[nmetro, nx] e_x_latent;
  vector[nmetro] e_index_latent;
}
model {
  {
    vector[nmetro] index_latent;
    matrix[nmetro, nx] x_latent;
    for(i in 1:nx){
      x_latent[,i] = mu_cs[i] + sigma_cs[i] * e_x_latent[,i];
    }
    index_latent = alpha_cs + x_latent * beta_cs + tau_cs * e_index_latent;

    to_vector(x_cs) ~ normal(to_vector(x_latent), to_vector(x_cs_se));
    index_cs ~ normal(index_latent, index_cs_se);
    e_index_latent ~ normal(0, 1);
    to_vector(e_x_latent) ~ normal(0, 1);

    alpha_cs ~ normal(alpha_cs_prior_mean, alpha_cs_prior_sd);
    beta_cs ~ normal(beta_cs_prior_mean, beta_cs_prior_sd);
    mu_cs ~ normal(mu_cs_prior_mean, mu_cs_prior_sd);
    tau_cs ~ normal(tau_cs_prior_mean, tau_cs_prior_sd);
    sigma_cs ~ normal(sigma_cs_prior_mean, sigma_cs_prior_sd);
  }
}
generated quantities {
  real alpha;
  vector[nx] beta;
  vector[nx] mu;
  vector<lower = 0>[nx] sigma;
  real<lower = 0> tau;

  beta = (beta_cs ./ sd_x) * sd_index;
  mu = mu_cs .* sd_x + mean_x;
  sigma = sd_x .* sigma_cs;
  tau = sd_index * tau_cs;
  alpha = sd_index * alpha_cs + mean_index - dot_product(mean_x, beta);
}
