data {
  int <lower=0> N;                // Number of samples
  int <lower=0> K;                // Number of taxa
  int <lower=0, upper=1> y[K, N]; // Binary abundance matrix
  vector[N] x;                    // Group
  vector[N] reads;                // Sequencing depth (log10 - centered)
}

parameters {
  vector[K] alpha;                // Intercepts for logistic regressions
  vector[K] beta_r;               // Parameter for reads
  vector[K] z;                    // Unscaled betas for logr
  real <lower=0, upper=1> nu;     // Assymmetry parametrer of Laplace prior
  real <lower=0> tau;             // Scale of Laplace prior
}

transformed parameters {
  vector[K] beta;                // DA estimates
  beta = z .* tau;
}

model {
  alpha ~ normal(0, 5);
  beta_r ~ normal(2, 2);
  nu ~ double_exponential(0.5, .05);
  tau ~ normal(0, 1);
  z ~ skew_double_exponential(0, 1, nu);

  for (k in 1:K) {
    y[k] ~ bernoulli_logit(alpha[k] + x * beta[k] + reads * beta_r[k]);
  }
}
