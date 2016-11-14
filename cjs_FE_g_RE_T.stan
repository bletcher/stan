// This models is derived from section 11.3 of "Stan Modeling Language
// User's Guide and Reference Manual"

functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }
  
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k;
      k = size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }
  
  matrix prob_uncaptured(int nInd, int nOcc,
                         matrix p, matrix phi) {
    matrix[nInd, nOcc] chi;
    
    for (i in 1:nInd) {
      chi[i, nOcc] = 1.0;
      for (t in 1:(nOcc - 1)) {
        int t_curr;
        int t_next;
        
        t_curr = nOcc - t;
        t_next = t_curr + 1;
        chi[i, t_curr] = (1 - phi[i, t_curr]) +
          phi[i, t_curr] *
          (1 - p[i, t_next - 1]) *
          chi[i, t_next];
      }
    }
    return chi;
  }
}

data {
  int<lower=0> nInd;            // Number of individuals
  int<lower=2> nOcc;            // Number of capture occasions
  int<lower=0,upper=1> y[nInd, nOcc];    // Capture-history

  int<lower=1> g;               // Number of groups
  int<lower=1,upper=g> riverOBear[nInd];     // Groups, 1 if in OBear, 0 otherwise
  real df;                      // Degree of freedom
  matrix[g,g] R;                // Scale matrix
}

transformed data {
  int<lower=0,upper=nOcc> first[nInd];
  int<lower=0,upper=nOcc> last[nInd];
  
  for (i in 1:nInd)
    first[i] = first_capture(y[i]);
  for (i in 1:nInd)
    last[i] = last_capture(y[i]);
}

parameters {
  vector<lower=0,upper=1>[g] mean_phi;   // Mean group-spec. survival
  real<lower=0,upper=1> p_g[g];          // Group-spec. recapture
  matrix[nOcc-1, g] eta_phi;
  cov_matrix[g] sigma;
}

transformed parameters {
  matrix<lower=0,upper=1>[nInd, nOcc - 1] phi;
  matrix<lower=0,upper=1>[nInd, nOcc - 1] p;
  matrix<lower=0,upper=1>[nInd, nOcc] chi;
  vector[g] mu_phi;
  
  // Constraints
//  mu = logit(mean_phi);

  for (u in 1:g)
    mu_phi[u] = logit(mean_phi[u]);

  for (i in 1:nInd) {
    for (t in 1:(first[i] - 1)) {
      phi[i, t] = 0;
      p[i, t] = 0;
    }
    for (t in first[i]:(nOcc - 1)) {
      phi[i, t] = inv_logit(eta_phi[t, riverOBear[i]]);
      p[i, t] = p_g[riverOBear[i]];
    }
  }
  
  chi = prob_uncaptured(nInd, nOcc, p, phi);
}

model {
  // Priors
  sigma ~ inv_wishart(df, R);

  for (t in 1:(nOcc - 1))
    increment_log_prob(multi_normal_log(row(eta_phi, t), mu_phi, sigma));
  mean_phi ~ uniform(0, 1);

  // for recapture parameters
  p_g ~ uniform(0, 1);
  
  // Likelihood
  for (i in 1:nInd) {
    if (first[i] > 0) {
      for (t in (first[i] + 1):last[i]) {
        1 ~ bernoulli(phi[i, t - 1]);
        y[i, t] ~ bernoulli(p[i, t - 1]);
      }
      1 ~ bernoulli(chi[i, last[i]]);
    }
  }
}

generated quantities {
//  real<lower=0> sigma2;
  
//  sigma2 = square(sigma);
}
