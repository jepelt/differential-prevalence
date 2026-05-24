data {
  int <lower=0> N;                // Number of samples
  int <lower=0> K;                // Number of taxa
  real y[K, N];                   // Relative abundance matrix
  vector[N] x;                    // Group
}

parameters {
  vector[K] alpha;                // Intercepts for logistic regressions
  vector[K] z;                    // Unscaled betas for logr
  vector <lower=0> [K] sigma;     // SDs of errors
  real <lower=0, upper=1> nu;     // Assymmetry parametrer of Laplace prior
  real <lower=0> tau;             // Scale of Laplace prior
}

transformed parameters {
  vector[K] beta;                 // DA estimates
  beta = z * tau;
}

model {
  alpha ~ normal(0, 5);
  sigma ~ gamma(1, 1);
  nu ~ double_exponential(0.5, .05);
  tau ~ normal(0, 1);
  z ~ skew_double_exponential(0, 1, nu);

  for (k in 1:K) {
    to_vector(y[k]) ~ normal(alpha[k] + x * beta[k], sigma[k]);
  }
}
