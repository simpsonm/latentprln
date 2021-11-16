data {
  int<lower = 1> nmetro;
  int<lower = 1> nx;
  vector[nmetro] index;
  vector<lower = 0>[nmetro] index_se;
  matrix[nmetro, nx] x;
}
transformed data {
  vector<lower = 0>[nmetro] index_var;
  index_var = index_se .* index_se;
}
parameters {
  real alpha;
  vector[nx] beta;
  real<lower = 0> sigma2;
}
model {
  index ~ normal(alpha + x*beta, sqrt(index_var + sigma2));
}
