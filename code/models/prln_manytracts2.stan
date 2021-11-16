data {
  int<lower = 2> nknot;
  ordered[nknot] knots;
  int<lower = 2> nbin;
  int<lower = 1, upper = nknot> pknot_idx_lb[nbin];
  int<lower = 1, upper = nknot> pknot_idx_ub[nbin];
  ordered[nbin - 1] bounds;
  int<lower = 1> ntract;
  matrix[ntract, nbin] bin_est;
  matrix<lower = 0>[ntract, nbin] bin_se;
  vector[ntract] mean_est;
  vector<lower = 0>[ntract] mean_se;
  vector[ntract] median_est;
  vector<lower = 0>[ntract] median_se;

  simplex[nknot] pknot_tract_prior_loc[ntract];
  real<lower = 0> pknot_tract_prior_scale;
  real alpha_prior_mean;
  real<lower = 0> alpha_prior_sd;  
}
transformed data {
  int median2knot[ntract]; // given a median est, gives the nearest knot from above
  vector[nknot - 1] knot_means; // note: last bin is always pareto
  vector[ntract] median_cdf_est;
  vector<lower = 0>[ntract] median_cdf_se;
  matrix[ntract, nbin + 2] tract_full_ests;
  matrix<lower = 0>[ntract, nbin + 2] tract_full_ses;
  int<lower = 1, upper = nknot - 1> nalpha[ntract];
  int<lower = 1, upper = nknot - 1> nunif[ntract];
  int totalalpha;

  for(i in 1:ntract){
    for(j in 1:nknot){
      if(knots[j] <= median_est[i]){
	median2knot[i] = j + 1;
      }
    }
    nunif[i] = median2knot[i] - 1;
    nalpha[i] = nknot - nunif[i];
  }
  totalalpha = sum(nalpha);

  for(i in 1:(nknot - 1)){
    knot_means[i] = (knots[i+1] + knots[i])/2;
  }

  median_cdf_est = rep_vector(0.5, ntract);
  tract_full_ests = append_col(append_col(bin_est, mean_est), median_cdf_est);

  for(i in 1:ntract){
    int bin;
    real nknotINbin;
    for(j in 1:nbin){
      if(pknot_idx_lb[j] <= nunif[i] && pknot_idx_ub[j] >= nunif[i]){
	  bin = j;
	  nknotINbin = pknot_idx_ub[j] - pknot_idx_lb[j] + 1.0;
      }
    }
    median_cdf_se[i] = (median_se[i] * bin_est[i, bin] / nknotINbin) /
      (knots[median2knot[i]] - knots[nunif[i]]);
  }
  tract_full_ses = append_col(append_col(bin_se, mean_se), median_cdf_se);
  
}
parameters {
  simplex[nknot] pknot_tract[ntract];
  vector<lower = 1>[totalalpha] alpha_tract;
}
transformed parameters {
  matrix[ntract, nbin + 2] tract_true_values;
  vector[totalalpha] upper_knot_means;  

  for(i in 1:ntract){
    for(j in 1:nbin){
      tract_true_values[i,j] = sum(pknot_tract[i, pknot_idx_lb[j]:pknot_idx_ub[j]]);
    }
  }

  {
    int idx;
    idx = 1;
    for(i in 1:ntract){
      for(j in 1:(nalpha[i] - 1)){
	real lb;
	real ub;
	real alpha;
	lb = knots[nunif[i] + j];
	ub = knots[nunif[i] + j + 1];
	alpha = alpha_tract[idx];
	upper_knot_means[idx] = (alpha / (alpha - 1)) * lb * ub *
	  (pow(ub, alpha - 1) - pow(lb, alpha - 1)) /
	  (pow(ub, alpha)     - pow(lb, alpha));
	idx = idx + 1;
      }
      upper_knot_means[idx] = knots[nknot] *
	(alpha_tract[idx] / (alpha_tract[idx] - 1));
      idx = idx + 1;
    }
  }

  {
    int idx1;
    int idx2;
    idx1 = 0;
    idx2 = 0;
    for(i in 1:ntract){
      idx2 = idx2 + nalpha[i];
      idx1 = idx2 - nalpha[i] + 1;
      // mean
      tract_true_values[i, nbin + 1] =
	dot_product(pknot_tract[i, 1:nunif[i]], knot_means[1:nunif[i]]) +
	dot_product(pknot_tract[i, (nunif[i] + 1):nknot], upper_knot_means[idx1:idx2]);
      // cdf at median
      tract_true_values[i, nbin + 2] = sum(pknot_tract[i, 1:nunif[i]]) -
	pknot_tract[i, nunif[i]] * (knots[median2knot[i]] - median_est[i]) /
	(knots[median2knot[i]] - knots[nunif[i]]);
      // note: cdf at median = sum of all uniform bin probs -
      //   bin prob [nunif] * (upper knot - median est) / (upper knot - lower knot)
    }
  }
}
model {
  to_vector(tract_full_ests) ~ normal(to_vector(tract_true_values),
				      to_vector(tract_full_ses));

  for(i in 1:ntract){
    pknot_tract[i] ~ dirichlet(pknot_tract_prior_loc[i] / pknot_tract_prior_scale);
  }
 

  alpha_tract ~ normal(alpha_prior_mean, alpha_prior_sd);
}
