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
  vector[nbin - 1] lower_knot_means; // note: last bin is always pareto
  vector[nbin + 1] tract_full_ests;
  vector<lower = 0>[nbin + 1] tract_full_ses;
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

  for(i in 1:(nbin - 1)){
    lower_knot_means[i] = (knots[i+1] + knots[i]) / 2.0;
  }

  tract_full_ests = append_row(bin_est, mean_est);
  tract_full_ses = append_row(bin_se, mean_se);

}
parameters {
  simplex[nbin] pknot;
  vector<lower = 1>[nalpha] alpha;
}
transformed parameters {
  vector[nest] tract_true_values;

  {
    // these are all intermediate quantities for computing
    // tract_true_values, and not intrinsically interesting
    // also, some elements of these are constant during 
    // sampling, which confuses Stan's warnings
    vector[nalpha] upper_knot_means;
    vector[nbin] pknot_times_muknot;

    // bin estimates
    tract_true_values[1:nbin] = pknot;


    // means of upper bins
    for(j in 1:(nalpha - 1)){
      real lb;
      real ub;
      lb = knots[nunif + j];
      ub = knots[nunif + j + 1];
      upper_knot_means[j] = (alpha[j] / (alpha[j] - 1)) * lb * ub *
	(pow(ub, alpha[j] - 1) - pow(lb, alpha[j] - 1)) /
	(pow(ub, alpha[j])     - pow(lb, alpha[j]));
    }
    upper_knot_means[nalpha] = knots[nbin] *
      (alpha[nalpha] / (alpha[nalpha] - 1));

    // pknot_times_muknot
    pknot_times_muknot[1:nunif] =
      pknot[1:nunif] .* lower_knot_means[1:nunif];
    pknot_times_muknot[(nunif + 1):nbin] =
      pknot[(nunif + 1):nbin] .* upper_knot_means[1:nalpha];

    // mean
    tract_true_values[nbin + 1] = sum(pknot_times_muknot);
  }
}
model {
  tract_full_ests ~ normal(tract_true_values, tract_full_ses);
  pknot ~ dirichlet(pknot_prior_loc / pknot_prior_scale);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
}
