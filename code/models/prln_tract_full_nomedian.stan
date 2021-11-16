data {
  int<lower = 2> nbin;
  ordered[nbin] knots;  
  int<lower = 1> nquant;    
  ordered[nquant] quantiles;
  int<lower = 1> nshare;
  vector<lower = 0>[nshare] share_lb;  
  vector[nshare] share_ub;  // < 0 implies share_ub = Infinity
  vector[nbin] bin_est;
  vector<lower = 0>[nbin] bin_se;
  vector[nquant] quant_est;
  vector<lower = 0>[nquant] quant_se;
  vector[nshare] share_est;
  vector<lower = 0>[nshare] share_se;  
  real<lower = 0, upper = 1> gini_est;
  real<lower = 0> gini_se;
  real mean_est;
  real<lower = 0> mean_se;
  simplex[nbin] pknot_prior_loc;
  real<lower = 0> pknot_prior_scale;
  real<lower = 1> alpha_prior_mean;
  real<lower = 0> alpha_prior_sd;  
}
transformed data {
  int nest; // total number of estimates
  int quantile2knot[nquant]; // nearest knot from below
  int sharelb2knot[nshare];  // nearest knot from below
  int shareub2knot[nshare];  // nearest knot from below
  int median2knot;           // nearest knot from above
  vector[nbin - 1] lower_knot_means; // note: last bin is always pareto
  vector[nbin - 1] lower_integrated_lorenz;
  vector[nquant] quant_cdf_est;
  vector<lower = 0>[nquant] quant_cdf_se;
  vector[nbin + 1 + nquant + nshare + 1] tract_full_ests;
  vector<lower = 0>[nbin + 1 + nquant + nshare + 1] tract_full_ses;
  int<lower = 1, upper = nbin - 1> nalpha;
  int<lower = 1, upper = nbin - 1> nunif;

  nest = nbin + 1 + nquant + nshare + 1;
  
  for(j in 1:nbin){
    if(knots[j] <= mean_est){
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

  for(j in 1:nbin){
    for(k in 1:nquant){
      if(knots[j] <= quant_est[k]){
	quantile2knot[k] = j;
      }
    }
  }

  for(j in 1:nbin){
    for(k in 1:nshare){
      if(knots[j] <= share_lb[k]){
	sharelb2knot[k] = j;
      }
      if(knots[j] <= share_ub[k]){
	shareub2knot[k] = j;
      }
    }
  }

  // if a share_ub < 0, that means that the real ub is infinity
  for(k in 1:nshare){
    if(share_ub[k] < 0.0){
      shareub2knot[k] = nbin;
    }
  }
  
  for(i in 1:(nbin - 1)){
    lower_knot_means[i] = (knots[i+1] + knots[i]) / 2.0;
    lower_integrated_lorenz[i] = (1 + knots[i]/(knots[i] + knots[i+1])) / 3.0;
  }

  for(j in 1:nquant){
    quant_cdf_est[j] = quantiles[j];
  }

  tract_full_ests =
    append_row(append_row(append_row(append_row(bin_est, mean_est),
				     quant_cdf_est),
			  share_est),
	       gini_est);

  for(k in 1:nquant){
    if(quantile2knot[k] <= nunif){
      quant_cdf_se[k] = (quant_se[k] * bin_est[quantile2knot[k]]) /
	(knots[quantile2knot[k] + 1] - knots[quantile2knot[k]]);
    } else if(quantile2knot[k] < nbin){
      real plower;
      real pupper;
      real lb;
      real ub;
      real myalpha;
      plower = sum(bin_est[quantile2knot[k]:nbin]) + 0.0001;
      pupper = plower - bin_est[quantile2knot[k]] + 0.0001;
      lb = knots[quantile2knot[k]];
      ub = knots[quantile2knot[k]+1];
      myalpha = log(plower/pupper)/log(ub/lb);
      
      if(myalpha <= 1){
	quant_cdf_se[k] = (quant_se[k] * bin_est[quantile2knot[k]]) /
	  (knots[quantile2knot[k] + 1] - knots[quantile2knot[k]]);
      } else {
	quant_cdf_se[k] = (quant_se[k] * bin_est[quantile2knot[k]]) *
	  ((myalpha / quant_est[k]) * pow(lb/quant_est[k], myalpha)) /
	  (1 - pow(lb/ub, myalpha));
      }
    } else {
      int idx;
      real plower;
      real pupper;
      real lb;
      real ub;
      real myalpha;
      idx = nbin - 1;
      myalpha = 0;
      while(myalpha <= 1 && idx > nunif){
	plower = sum(bin_est[idx:nbin]) + 0.0001;
	pupper = plower - bin_est[idx] + 0.0001;
	lb = knots[idx];
	ub = knots[idx+1];
	myalpha = log(plower/pupper)/log(ub/lb);
	idx = idx - 1;
      }
      if(myalpha <= 1){
	myalpha = 1.0001;
      }
      quant_cdf_se[k] = (quant_se[k] * bin_est[nbin]) *
	(myalpha / quant_est[k]) * pow(knots[nbin]/quant_est[k], myalpha);
    }
  }

  tract_full_ses =
    append_row(append_row(append_row(append_row(bin_se, mean_se),
				     quant_cdf_se),
			  share_se),
	       gini_se);
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
    vector[nshare] lorenzk_lb; 
    vector[nshare] lorenzk_ub;

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

    // pknot_times_muknot, for gini and shares
    pknot_times_muknot[1:nunif] =
      pknot[1:nunif] .* lower_knot_means[1:nunif];
    pknot_times_muknot[(nunif + 1):nbin] =
      pknot[(nunif + 1):nbin] .* upper_knot_means[1:nalpha];

    // within bin lorenz curves
    for(j in 1:nshare){
      // lower bound
      int idx = sharelb2knot[j];
      if(idx <= nunif){
	lorenzk_lb[j] = (share_lb[j] * share_lb[j] - knots[idx] * knots[idx]) /
	  ((knots[idx + 1] - knots[idx]) * 2 * lower_knot_means[idx]);
      } else if(idx < nbin){
	real myalpha = alpha[idx - nunif];
	lorenzk_lb[j] = (myalpha / (myalpha - 1)) *
	  (knots[idx] / upper_knot_means[idx - nunif]) *
	  (1 - pow(knots[idx]/share_lb[j], myalpha - 1)) /
	  (1 - pow(knots[idx]/knots[idx + 1], myalpha));
      } else {
	real myalpha = alpha[nalpha];
	lorenzk_lb[j] = (myalpha / ((myalpha - 1) * upper_knot_means[nalpha])) *
	  (knots[nbin] - pow(knots[nbin], myalpha) * pow(share_lb[j], 1 - myalpha));
      }
      // upper bound
      idx = shareub2knot[j];
      if(idx <= nunif){
	lorenzk_ub[j] = (share_ub[j] * share_ub[j] - knots[idx] * knots[idx]) /
	  ((knots[idx + 1] - knots[idx]) * 2 * lower_knot_means[idx]);
      } else if(idx < nbin){
	real myalpha = alpha[idx - nunif];
	lorenzk_ub[j] = (myalpha / (myalpha - 1)) *
	  (knots[idx] / upper_knot_means[idx - nunif]) *
	  (1 - pow(knots[idx]/share_ub[j], myalpha - 1)) /
	  (1 - pow(knots[idx]/knots[idx + 1], myalpha));
      } else {
	if(share_ub[j] > 0){
	  real myalpha = alpha[nalpha];
	  lorenzk_ub[j] =
	    (myalpha / ((myalpha - 1) * upper_knot_means[nalpha])) *
	    (knots[nbin] - pow(knots[nbin], myalpha) *
	     pow(share_ub[j], 1 - myalpha));
	} else {
	  lorenzk_ub[j] = 1.0;
	}
      }
    }

    // pmui, for Gini
    pmui[1:nunif] = pknot_times_muknot[1:nunif] .* 
      (pknot[1:nunif] .* lower_integrated_lorenz[1:nunif] +
       sum_prob_above[1:nunif]);
    pmui[(nunif + 1):nbin] = pknot_times_muknot[(nunif + 1):nbin] .* 
      (pknot[(nunif + 1):nbin] .* upper_integrated_lorenz[1:nalpha] +
       sum_prob_above[(nunif + 1):nbin]);

    // mean
    tract_true_values[nbin + 1] = sum(pknot_times_muknot);

    // cdf at quantile
    for(j in 1:nquant){
      int idx1;
      int idx2;
      idx1 = quantile2knot[j];
      idx2 = idx1 + 1;
      if(idx1 <= nunif){
	tract_true_values[nbin + 1 + j] = sum(pknot[1:idx1]) -
	  pknot[idx1] * (knots[idx2] - quant_est[j]) /
	  (knots[idx2] - knots[idx1]);
	// note: cdf at quant = sum of bin probs up to and including this bin -
	//   bin prob [this bin] * (upper knot - quant est) / (upper knot - lower knot)
      } else if(idx1 < nbin){
	real myalpha = alpha[idx1 - nunif];
	tract_true_values[nbin + 1 + j] = sum(pknot[1:(idx1-1)]) +
	  pknot[idx1] * (1 - pow(knots[idx1]/quant_est[j], myalpha)) /
	  (1 - pow(knots[idx1]/knots[idx2], myalpha));
      } else {
	real myalpha = alpha[nalpha];
	tract_true_values[nbin + 1 + j] = sum(pknot[1:(idx1-1)]) +
	  pknot[idx1] * (1 - pow(knots[idx1]/quant_est[j], myalpha));	
      }
    }

    // income shares
    for(j in 1:nshare){
      tract_true_values[nbin + 1 + nquant + j] =
	(sum(pknot_times_muknot[sharelb2knot[j]:shareub2knot[j]]) -
	 pknot_times_muknot[shareub2knot[j]] * (1 - lorenzk_ub[j]) -
	 pknot_times_muknot[sharelb2knot[j]] * lorenzk_lb[j]) /
	tract_true_values[nbin + 1];
    }

    // gini coefficient
    tract_true_values[nbin + 1 + nquant + nshare + 1] =
      1 - 2*sum(pmui)/tract_true_values[nbin + 1];
  }
}
model {
  tract_full_ests ~ normal(tract_true_values, tract_full_ses);
  pknot ~ dirichlet(pknot_prior_loc / pknot_prior_scale);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
}
