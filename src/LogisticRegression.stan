data {
  int<lower=1> N; //number of responses
  int<lower=0> d; //number of fixed effect covariates (not including intercept)
  matrix[N,d] X; //Fixed efffects design matrix.
  int<lower=0,upper=1> y[N]; //response
  real<lower=0> sigmab; //prior variance of beta (fixed effects)
  real mualpha; //prior mean of alpha (intercept)
  real<lower=0> sigmaalpha; //prior sd of alpha (intecept)
}

parameters {
  real alpha;  //logistic regression intercept 
  vector[d] beta; //logistic regression coefficients for fixed effect
}

model {
  vector[N] mu_fixed;
  if (d > 0) {
    mu_fixed = alpha + X * beta;
  } else {
    mu_fixed = rep_vector(alpha, N);    
  }
  y ~ bernoulli_logit(mu_fixed);
  for (i in 1:d) {
    beta[i] ~ normal(0, sigmab);
  }
  alpha ~ normal(mualpha, sigmaalpha); 
}