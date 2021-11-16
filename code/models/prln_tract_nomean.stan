data {
  int<lower = 2> nbin;
  ordered[nbin] knots;  
  vector[nbin] bin_est;
  vector<lower = 0>[nbin] bin_se;
  real mean_est;
  real<lower = 0> mean_se;
  simplex[nbin] pknot_prior_loc;
  real<lower = 0> pknot_prior_scale;
  real<lower = 1> alpha_prior_mean;
  real<lower = 0> alpha_prior_sd;  
}
transformed data {
  int nest; // total number of estimates
  // given a median estimate, gives the nearest knot from above
  int median2knot; 
  real median_cdf_est;
  real<lower = 0> median_cdf_se;
  vector[nbin + 1] tract_full_ests;
  vector<lower = 0>[nbin + 1] tract_full_ses;
  int<lower = 1, upper = nbin - 1> nalpha;
  int<lower = 1, upper = nbin - 1> nunif;

  nest = nbin + 2;
  
  for(j in 1:nbin){
    if(knots[j] <= median_est){
      median2knot = j + 1;
    }
  }
  nunif = median2knot - 1;

  // last bin must always be pareto
  if(nunif == nbin){
    nunif = nunif - 1;
    median2knot = median2knot - 1;
  }
  
  nalpha = nbin - nunif;

  median_cdf_est = 0.5;
  tract_full_ests = append_row(bin_est, median_cdf_est);

  median_cdf_se = median_se * bin_est[nunif] / (knots[nunif + 1] - knots[nunif]);
  tract_full_ses = append_row(bin_se, median_cdf_se);
}
parameters {
  simplex[nbin] pknot;
  vector<lower = 1>[nalpha] alpha;
}
transformed parameters {
  vector[nest] tract_true_values;

  // bin estimates
  tract_true_values[1:nbin] = pknot;

  // cdf at median
  // note: cdf at median = sum of bin probs up to and including this bin minus
  //   bin prob [this bin] * (upper knot - median est) / (upper knot - lower knot)
  tract_true_values[nbin + 1] =
    sum(pknot[1:nunif]) - pknot[nunif] *
    (knots[median2knot] - median_est) /
    (knots[median2knot] - knots[nunif]);
}
model {
  tract_full_ests ~ normal(tract_true_values, tract_full_ses);
  pknot ~ dirichlet(pknot_prior_loc / pknot_prior_scale);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
}
