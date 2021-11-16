data {
  int<lower = 2> nbin;
  ordered[nbin] knots;  
  vector[nbin] bin_est;
  vector<lower = 0>[nbin] bin_se;
  real<lower = 0, upper = 1> gini_est;
  real<lower = 0> gini_se;
  real mean_est;
  real<lower = 0> mean_se;
  real median_est;
  real<lower = 0> median_se;  
  simplex[nbin] pknot_prior_loc;
  real<lower = 0> pknot_prior_scale;
  real<lower = 1> alpha_prior_mean;
  real<lower = 0> alpha_prior_sd;  
}
transformed data {
  int nest; // total number of estimates
  // given a median estimate, gives the nearest knot from above
  int median2knot; 
  vector[nbin - 1] lower_knot_means; // note: last bin is always pareto
  vector[nbin - 1] lower_integrated_lorenz;
  real median_cdf_est;
  real<lower = 0> median_cdf_se;
  vector[nbin + 2 + 1] tract_full_ests;
  vector<lower = 0>[nbin + 2 + 1] tract_full_ses;
  int<lower = 1, upper = nbin - 1> nalpha;
  int<lower = 1, upper = nbin - 1> nunif;

  nest = nbin + 2 + 1;
  
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

  for(i in 1:(nbin - 1)){
    lower_knot_means[i] = (knots[i+1] + knots[i]) / 2.0;
    lower_integrated_lorenz[i] = (1 + knots[i]/(knots[i] + knots[i+1])) / 3.0;
  }

  median_cdf_est = 0.5;
  tract_full_ests = append_row(append_row(append_row(bin_est, mean_est),
					  median_cdf_est), gini_est);

  median_cdf_se = median_se * bin_est[nunif] / (knots[nunif + 1] - knots[nunif]);

  tract_full_ses = append_row(append_row(append_row(bin_se, mean_se),
					 median_cdf_se), gini_se);
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
    vector[nalpha] upper_integrated_lorenz;
    vector[nbin] pknot_times_muknot;
    vector[nbin] pmui; // sum(pmui) / mu = integrated lorenz curve for the tract
    vector[nbin] sum_prob_above;

    // bin estimates
    tract_true_values[1:nbin] = pknot;

    // sum of all probabilities above current bin, for Gini
    sum_prob_above[nbin] = 0.0;
    for(j in 1:(nbin - 1)){
      sum_prob_above[nbin - j] = sum_prob_above[nbin - j + 1] + pknot[nbin - j + 1];
    }

    // means of upper bins
    for(j in 1:(nalpha - 1)){
      real lb;
      real ub;
      lb = knots[nunif + j];
      ub = knots[nunif + j + 1];
      upper_knot_means[j] = (alpha[j] / (alpha[j] - 1)) * lb * ub *
	(pow(ub, alpha[j] - 1) - pow(lb, alpha[j] - 1)) /
	(pow(ub, alpha[j])     - pow(lb, alpha[j]));
      upper_integrated_lorenz[j] =
	1 - ((alpha[j]*pow(ub, 1 - alpha[j])) / (2*alpha[j] - 1)) *
	(pow(ub, 2*alpha[j] - 1) - pow(lb, 2*alpha[j] - 1)) /
	(pow(ub, alpha[j]) - pow(lb, alpha[j]));
    }
    upper_knot_means[nalpha] = knots[nbin] *
      (alpha[nalpha] / (alpha[nalpha] - 1));
    upper_integrated_lorenz[nalpha] = 1 - alpha[nalpha] / (2*alpha[nalpha] - 1);

    // pknot_times_muknot, for Gini and mean
    pknot_times_muknot[1:nunif] =
      pknot[1:nunif] .* lower_knot_means[1:nunif];
    pknot_times_muknot[(nunif + 1):nbin] =
      pknot[(nunif + 1):nbin] .* upper_knot_means[1:nalpha];

    // pmui, for Gini
    pmui[1:nunif] = pknot_times_muknot[1:nunif] .* 
      (pknot[1:nunif] .* lower_integrated_lorenz[1:nunif] +
       sum_prob_above[1:nunif]);
    pmui[(nunif + 1):nbin] = pknot_times_muknot[(nunif + 1):nbin] .* 
      (pknot[(nunif + 1):nbin] .* upper_integrated_lorenz[1:nalpha] +
       sum_prob_above[(nunif + 1):nbin]);

    // mean
    tract_true_values[nbin + 1] = sum(pknot_times_muknot);

    // cdf at median
    // note: cdf at median = sum of bin probs up to and including this bin -
    //   bin prob [this bin] * (upper knot - median est) / (upper knot - lower knot)
    tract_true_values[nbin + 2] =
      sum(pknot[1:nunif]) - pknot[nunif] *
      (knots[median2knot] - median_est) /
      (knots[median2knot] - knots[nunif]);
  
    // gini coefficient
    tract_true_values[nbin + 2 + 1] =
      1 - 2*sum(pmui)/tract_true_values[nbin + 1];
  }
}
model {
  tract_full_ests ~ normal(tract_true_values, tract_full_ses);
  pknot ~ dirichlet(pknot_prior_loc / pknot_prior_scale);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
}
