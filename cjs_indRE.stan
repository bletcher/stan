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
  real<lower=0,upper=1> mean_phi;       // Mean survival
  real<lower=0,upper=1> mean_p;         // Mean recapture 
  real epsilon[nInd];
  real<lower=0> sigma;
}

transformed parameters {
  matrix<lower=0,upper=1>[nInd, nOcc - 1] phi;
  matrix<lower=0,upper=1>[nInd, nOcc - 1] p;
  matrix<lower=0,upper=1>[nInd, nOcc] chi;
  real mu;
  
  // Constraints
  mu = logit(mean_phi);
  for (i in 1:nInd) {
    for (t in 1:(first[i] - 1)) {
      phi[i, t] = 0;
      p[i, t] = 0;
    }
    for (t in first[i]:(nOcc - 1)) {
      phi[i, t] = inv_logit(mu + epsilon[i]);
      p[i, t] = mean_p;
    }
  }
  
  chi = prob_uncaptured(nInd, nOcc, p, phi);
}

model {
  // Priors
  epsilon ~ normal(0, sigma);
  mean_phi ~ uniform(0, 1);
  sigma ~ uniform(0, 1); //cauchy(0,1)
  mean_p ~ uniform(0, 1);
  
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
  real<lower=0> sigma2;
  
  sigma2 = square(sigma);
}
