data {
  int<lower = 2> nbin;
  ordered[nbin] knots;  
  vector[nbin] bin_est;
  vector<lower = 0>[nbin] bin_se;
  simplex[nbin] pknot_prior_loc;
  real<lower = 0> pknot_prior_scale;
  real<lower = 1> alpha_prior_mean;
  real<lower = 0> alpha_prior_sd;  
}
transformed data {
  int nest; // total number of estimates
  int<lower = 1, upper = nbin - 1> nalpha;
  int<lower = 1, upper = nbin - 1> nunif;

  nest = nbin + 1;

  {
    vector[nbin] binsum;
    binsum = cumulative_sum(bin_est);
    nunif = 1;
    for(j in nbin:1){
      if(binsum[j] < 0.5){
	nunif = j + 1;
      }
    }
  }

  // last bin must always be pareto
  if(nunif == nbin){
    nunif = nunif - 1;
  }
  
  nalpha = nbin - nunif;
}
parameters {
  simplex[nbin] pknot;
  vector<lower = 1>[nalpha] alpha;
}
model {
  bin_est ~ normal(pknot, bin_se);
  pknot ~ dirichlet(pknot_prior_loc / pknot_prior_scale);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
}
