data {
  int<lower = 1> nmetro;
  int<lower = 1> nx;
  vector[nmetro] index;
  vector<lower = 0>[nmetro] index_se;
  matrix[nmetro, nx] x;
  matrix<lower = 0>[nmetro, nx] x_se;
  real alpha_prior_mean;
  real<lower = 0> alpha_prior_sd;
  vector[nx] beta_prior_mean;
  vector<lower = 0>[nx] beta_prior_sd;
  vector[nmetro] mu_prior_mean;
  vector<lower = 0>[nmetro] mu_prior_sd;
  
}
parameters {
  real alpha;
  vector[nx] beta;
  vector[nx] mu;
  vector<lower = 0>[nx] sigma;
  real<lower = 0> tau;
  matrix[nmetro, nx] x_latent;
  vector[nmetro] index_latent;
}
model {
  to_vector(x) ~ normal(to_vector(x_latent), to_vector(x_se));
  index ~ normal(index_latent, index_se);
  index_latent ~ normal(alpha + x_latent * beta, tau);
  for(i in 1:nx){
    x_latent[,i] ~ normal(mu[i], sigma[i]);
  }
}
