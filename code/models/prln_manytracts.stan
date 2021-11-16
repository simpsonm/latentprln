// This version of the model allows for:
//  1) knots which are different from the bin boundaries, and
//  2) fitting multiple tracts at the same time
// Otherwise it is the same as the other models

// note: if a bin has bounds LB and UB, then the bin is [LB, UB)
//       so a knot is contained in a bin if LB <= knot < UB
// note: the first knot is the left knot of the first knot-bin
functions {
  real compute_bin_prob_nknot0(int nknot, int nunif, int nbin,
			       int idx_lb, int idx_ub, int bin_id,
			       data vector bounds, data vector knots,
			       vector pknot_tract, vector alpha_tract){
    real out;

    if(idx_lb <= nunif){
      // the knot-bin is uniform
      out = pknot_tract[idx_lb] *
	(bounds[bin_id + 1] - bounds[bin_id]) /
	(knots[idx_lb + 1] - knots[idx_lb]);
    } else if(idx_lb < nknot){
      // the knot-bin is truncated pareto
      real alpha;
      real knot_lb;
      real knot_ub;
      real pbounds;
      knot_lb = knots[idx_lb];
      knot_ub = knots[idx_lb + 1];
      alpha = alpha_tract[idx_lb - nunif];
      pbounds = pow(knot_lb / bounds[bin_id], alpha) -
	pow(knot_lb / bounds[bin_id + 1], alpha);
      out = pknot_tract[idx_lb] * pbounds / pareto_cdf(knot_ub, knot_lb, alpha);
    } else {
      // the knot-bin is untruncated pareto
      real alpha;
      real knot_lb;
      real pbounds;
      knot_lb = knots[idx_lb];
      alpha = alpha_tract[idx_lb - nunif];
      if(bin_id < nbin){
	// bin is not the rightmost, open ended bin
	pbounds = pow(knot_lb / bounds[bin_id], alpha) -
	  pow(knot_lb / bounds[bin_id + 1], alpha);
      } else {
	// bin is the rightmost, open ended bin
	pbounds = pow(knot_lb / bounds[bin_id], alpha);
      }
      out = pknot_tract[idx_lb] * pbounds;
    }
    return out;
  }

  // note: this will also compute the portion of the bin probability associated
  //       with the left knot-bin and right knot-bin if there are > 1 knots in the bin
  real compute_bin_prob_nknot1(int nknot, int nunif, int nbin,
			       int idx_lb, int idx_ub, int bin_id,
			       data vector bounds, data vector knots,
			       vector pknot_tract, vector alpha_tract){
    real left_prob;
    real right_prob;

    // first do the left knot-bin probability
    if(idx_lb <= nunif){
      // the knot-bin is uniform
      if(knots[idx_lb + 1] == bounds[bin_id]){
	//  the knot and the boundary coincide
	//   need to manually set this probability to zero so that we
	//       avoid infinite gradients generated by autodiff
	//       (which branch of this if statement we enter does NOT
	//        depend on unknown parameters, so it doesn't break HMC)	
	left_prob = 0;
      } else {
	left_prob = pknot_tract[idx_lb] *
	  (knots[idx_lb + 1] - bounds[bin_id]) /
	  (knots[idx_lb + 1] - knots[idx_lb]);
      }
    } else if(idx_lb < nknot){
      // the knot-bin is truncated pareto
      real alpha;
      real knot_lb;
      real knot_ub;
      knot_lb = knots[idx_lb];
      knot_ub = knots[idx_lb + 1];
      alpha = alpha_tract[idx_lb - nunif];
      if(knot_ub == bounds[bin_id]){
	//  the knot and the boundary coincide
	//   need to manually set this probability to zero so that we
	//       avoid infinite gradients generated by autodiff
	//       (which branch of this if statement we enter does NOT
	//        depend on unknown parameters, so it doesn't break HMC)	
	left_prob = 0;
      } else if(knot_lb == bounds[bin_id]) {
	//  the knot and the boundary coincide
	//   need to manually set this probability so that we
	//       avoid infinite gradients generated by autodiff
	//       (which branch of this if statement we enter does NOT
	//        depend on unknown parameters, so it doesn't break HMC)	
	left_prob = pknot_tract[idx_lb];
      } else {
	left_prob = pknot_tract[idx_lb] *
	  ( 1 - (pareto_cdf(bounds[bin_id], knot_lb, alpha) / 
		 pareto_cdf(knot_ub,        knot_lb, alpha)  ) );
      }
    } else {
      reject("compute_bin_prob_nknot1: idx_lb must be < nknot;",
    	     "found idx_lb = ", idx_lb,
    	     " bin idx = ", bin_id);
    }

    // now do the right knot-bin probability
    if(idx_ub - 1 <= nunif){
      // the knot-bin is uniform
      if(knots[idx_ub - 1] == bounds[bin_id + 1]){
	//  the knot and the boundary coincide
	//   need to manually set this probability to zero so that we
	//       avoid infinite gradients generated by autodiff
	//       (which branch of this if statement we enter does NOT
	//        depend on unknown parameters, so it doesn't break HMC)	
	right_prob = 0;
      } else {
	right_prob = pknot_tract[idx_ub - 1] *
	  (bounds[bin_id + 1] - knots[idx_ub - 1]) /
	  (knots[idx_ub] - knots[idx_ub - 1]);
      }
    } else if(idx_ub - 1 < nknot){
      // the knot-bin is truncated pareto
      real alpha;
      real knot_lb;
      real knot_ub;
      knot_lb = knots[idx_ub - 1];
      knot_ub = knots[idx_ub];
      alpha = alpha_tract[idx_ub - 1 - nunif];
      if(knot_lb == bounds[bin_id + 1]){
	//  the knot and the boundary coincide
	//   need to manually set this probability to zero so that we
	//       avoid infinite gradients generated by autodiff
	//       (which branch of this if statement we enter does NOT
	//        depend on unknown parameters, so it doesn't break HMC)	
	right_prob = 0;
      } else {
	// the knot is different from the boundary
	right_prob = pknot_tract[idx_ub - 1] *
	  pareto_cdf(bounds[bin_id + 1], knot_lb, alpha) /
	  pareto_cdf(knot_ub,            knot_lb, alpha);
      }
    } else {
      // the knot-bin is untruncated pareto (rightmost knot-bin)
      real alpha;
      real knot_lb;
      knot_lb = knots[idx_ub - 1];
      alpha = alpha_tract[idx_ub - 1 - nunif];
      if(bin_id == nbin){
	// we're in the last bin, so we get all of the last knot-bin
	right_prob = pknot_tract[idx_ub - 1];
      } else if(knot_lb == bounds[bin_id + 1]){
	// we're not in the last bin, but the knot and the boundary coincide
	//   need to manually set this probability to zero so that we
	//       avoid infinite gradients generated by autodiff
	//       (which branch of this if statement we enter does NOT
	//        depend on unknown parameters, so it doesn't break HMC)
	right_prob = 0;
      } else {
	// we're not in the last bin, and the knot is different from the boundary
	right_prob = pknot_tract[idx_ub - 1] *
	  pareto_cdf(bounds[bin_id + 1], knot_lb, alpha);
      }
    }

    return left_prob + right_prob;
  }

  real compute_bin_prob_nknotn(int nknot, int nunif, int nbin,
			       int idx_lb, int idx_ub, int bin_id,
			       data vector bounds, data vector knots,
			       vector pknot_tract, vector alpha_tract){
    real outside_prob;
    real middle_prob;
    outside_prob =
      compute_bin_prob_nknot1(nknot, nunif, nbin, idx_lb, idx_ub, bin_id,
			      bounds, knots, pknot_tract, alpha_tract);
    if(idx_lb + 1 > idx_ub - 2){
      reject("idx_lb + 1 > idx_ub - 2; ",
	     "idx_lb: ", idx_lb,
	     " idx_ub: ", idx_ub);
    }
    middle_prob = sum(pknot_tract[(idx_lb + 1):(idx_ub - 2)]);
    return outside_prob + middle_prob;
  }
}
data {
  int<lower = 1> ntract;
  int<lower = 2> nbin;
  int<lower = 2> nknot;                     // excluding the last knot (= infinity)
  int<lower = 2, upper = nknot - 2> nalpha[ntract]; // number of pareto knot-bins

  ordered[nbin]  bounds;                   // bounds of the bin estimates, w/ first = 0
  ordered[nknot] knots[ntract];            // knots for each tract
  matrix[ntract, nbin] bin_est;
  matrix<lower = 0>[ntract, nbin] bin_se;
  vector[ntract] mean_est;
  vector<lower = 0>[ntract] mean_se;
  vector[ntract] median_est;
  vector<lower = 0>[ntract] median_se;
  simplex[nknot] pknot_tract_prior_loc[ntract];
  real<lower = 0> pknot_tract_prior_scale;
  vector[sum(nalpha)] alpha_prior_mean;
  real<lower = 0> alpha_prior_sd;  
}
transformed data {
  int<lower = 2, upper = nknot - 2> nunif[ntract];        // number of uniform knot-bins
  vector[ntract*nknot - sum(nalpha)] lower_knot_means;    // means of the uniform bins
  int<lower = 1, upper = nbin> knot2bin[ntract, nknot];   // bin each knot is located in
  int<lower = 1, upper = nknot>     bin2knot_lb[ntract, nbin]; // max knot <= bin
  int<lower = 1, upper = nknot + 1> bin2knot_ub[ntract, nbin]; // min knot >= bin 
  int<lower = 0, upper = nknot> nknot_in_bin[ntract, nbin];  
  vector[ntract] median_cdf_est;
  vector<lower = 0>[ntract] median_cdf_se;
  matrix[ntract, nbin + 2] tract_full_ests;
  matrix<lower = 0>[ntract, nbin + 2] tract_full_ses;

  for(i in 1:ntract){
    nunif[i] = nknot - nalpha[i];
  }

  // find the index of the largest bound <= each knot
  // i.e. the index of the bin each knot is located in
  for(i in 1:ntract){
    for(j in 1:nknot){
      for(k in 1:nbin){
	if(knots[i,j] >= bounds[k]){
	  knot2bin[i,j] = k;  
	}
      }
    }
  }

  // find the index of the largest knot <= each bound
  // i.e. the index of the largest knot < or on the lower boundardy of each bin
  for(i in 1:ntract){
    for(j in 1:nbin){
      for(k in 1:nknot){
	if(bounds[j] >= knots[i,k]){
	  bin2knot_lb[i,j] = k;  
	}
      }
    }
  }

  // find the index of the smallest knot >= each bound
  // i.e. the index of hte smallest knot >= each bin
  // note: bin2knot_ub[i, j] = nknot + 1 indicates that the only knot above or
  //       or on the boundary of the bin is the last knot equal to infinity
  for(i in 1:ntract){
    bin2knot_ub[i, nbin] = nknot + 1;
    for(j in 1:(nbin - 1)){
      bin2knot_ub[i,j] = 1;
      for(k in 1:nknot){
	if(bounds[j + 1] > knots[i,k]){
	  bin2knot_ub[i,j] = bin2knot_ub[i,j] + 1;  
	}
      }
    }
  }

  // count the number of knots within each bin
  for(i in 1:ntract){
    for(j in 1:nbin){
      nknot_in_bin[i, j] = bin2knot_ub[i, j] - bin2knot_lb[i, j] - 1;
    }
  }

  // compute the means of the uniform knot-bins
  for(i in 1:ntract){
      for(j in 1:nunif[i]){
	int id;
	id = sum(nunif[1:i]) - nunif[i] + j;
	lower_knot_means[id] = (knots[i,j+1] + knots[i,j])/2;
      }
  }

  median_cdf_est = rep_vector(0.5, ntract);
  tract_full_ests = append_col(append_col(bin_est, mean_est), median_cdf_est);

  // uses the delta method and an approximation of the density as uniform within
  //   the bounds defining the bin containing the median, taking the associated
  //   bin estimate as the true probability within that bin
  // assumes the median estimate is within the last uniform knot-bin for each tract
  for(i in 1:ntract){
    int bin_id;
    bin_id = knot2bin[i,nunif[i]];
    if(bin_id < nbin){
      median_cdf_se[i] = median_se[i] * bin_est[i, bin_id] /
	(bounds[bin_id + 1] - bounds[bin_id]);
    } else {
      median_cdf_se[i] = median_se[i] * bin_est[i, bin_id] /
	(knots[i,nunif[i]] - bounds[bin_id]);
    }

  }
  tract_full_ses = append_col(append_col(bin_se, mean_se), median_cdf_se);
  
}
parameters {
  simplex[nknot] pknot_tract[ntract];
  vector<lower = 1>[sum(nalpha)] alpha_tract;
}
transformed parameters {
  matrix[ntract, nbin + 2] tract_true_values;
  vector[sum(nalpha)] upper_knot_means;

  for(i in 1:ntract){
    int alpha_id_lb;
    int alpha_id_ub;
    alpha_id_lb = sum(nalpha[1:i]) - nalpha[i] + 1;
    alpha_id_ub = sum(nalpha[1:i]);
    
    // means in the pareto bins, for constructing the tract means
    {
      real lb;
      real ub;
      real alpha;
      for(j in 1:(nalpha[i] - 1)){
	lb = knots[i, nunif[i] + j];
 	ub = knots[i, nunif[i] + j + 1];
	alpha = alpha_tract[sum(nalpha[1:i]) - nalpha[i] + j];
	upper_knot_means[sum(nalpha[1:i]) - nalpha[i] + j] =
	  (alpha / (alpha - 1)) * lb *
	  pareto_cdf(ub, lb, alpha - 1) / pareto_cdf(ub, lb, alpha);
      }
      lb = knots[i, nknot];
      alpha = alpha_tract[sum(nalpha[1:i])];
      upper_knot_means[sum(nalpha[1:i])] = (alpha / (alpha - 1)) * lb;
    }

    // bin probabilities
    {
      for(j in 1:nbin){
	int idx_lb;
	int idx_ub;
	idx_lb = bin2knot_lb[i,j];
	idx_ub = bin2knot_ub[i,j];
	if(0 == nknot_in_bin[i,j]){
	  tract_true_values[i, j] = 
	    compute_bin_prob_nknot0(nknot, nunif[i], nbin, idx_lb, idx_ub, j,
				    bounds, knots[i], pknot_tract[i],
				    alpha_tract[alpha_id_lb:alpha_id_ub]);
	} else if(1 == nknot_in_bin[i,j]){
	  tract_true_values[i, j] = 
	    compute_bin_prob_nknot1(nknot, nunif[i], nbin, idx_lb, idx_ub, j,
				    bounds, knots[i], pknot_tract[i],
				    alpha_tract[alpha_id_lb:alpha_id_ub]);
	} else {
	  tract_true_values[i, j] = 
	    compute_bin_prob_nknotn(nknot, nunif[i], nbin, idx_lb, idx_ub, j,
				    bounds, knots[i], pknot_tract[i],
				    alpha_tract[alpha_id_lb:alpha_id_ub]);

	}
      }
    }

    // mean
    tract_true_values[i, nbin + 1] =
      dot_product(pknot_tract[i, 1:nunif[i]],
		  lower_knot_means[(sum(nunif[1:i]) - nunif[i] + 1):sum(nunif[1:i])]) + 
      dot_product(pknot_tract[i, (nunif[i] + 1):nknot],
		  upper_knot_means[alpha_id_lb:alpha_id_ub]);

    // cdf at median
    tract_true_values[i, nbin + 2] = sum(pknot_tract[i, 1:(nunif[i] - 1)]) +
      pknot_tract[i, nunif[i]] * (median_est[i] - knots[i, nunif[i]]) /
      (knots[i, nunif[i] + 1] - knots[i, nunif[i]]);
    // note: cdf at median = sum of all uniform knot-bin probs except the last +
    //   bin prob [nunif] * (median est - lower knot) / (upper knot - lower knot)
    //   (we assume that the median is in the nunif knot-bin)
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
