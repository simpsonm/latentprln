data {
  int<lower = 1> nmetro;
  int<lower = 1> nx;
  vector[nmetro] index;
  vector<lower = 0>[nmetro] index_se;
  matrix[nmetro, nx] x;
  matrix<lower = 0>[nmetro, nx] x_se;
}
transformed data {
  matrix<lower = 0>[nmetro, nx] x_var;
  vector<lower = 0>[nmetro] index_var;
  x_var = x_se .* x_se;
  index_var = index_se .* index_se;
}
parameters {
  real alpha;
  vector[nx] beta;
  vector[nx] mu;
  vector<lower = 0>[nx] sigma2;
  real<lower = 0> tau2;
}
model {
  for(i in 1:nx){
    x[,i] ~ normal(mu[i], sqrt(sigma2[i] + x_var[,i]));
  }
  for(i in 1:nmetro){
    real mubar;
    real sig2bar;
    vector[nx] lambda;
    lambda = sigma2 ./ (sigma2 + to_vector(x_var[i]) );
    mubar = alpha + dot_product(mu + lambda .* (to_vector(x[i]) - mu), beta);
    sig2bar = index_var[i] + tau2 + dot_product(sigma2 .* (1 - lambda), beta .* beta);
    index ~ normal(mubar, sqrt(sig2bar));
  }
}
