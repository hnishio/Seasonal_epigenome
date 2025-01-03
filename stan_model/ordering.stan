
data {
  int<lower=0> N_constraints;
  int<lower=1, upper=10> var1[N_constraints];
  int<lower=1, upper=10> var2[N_constraints];
  real<lower=0> sigma; 
}

parameters {
  vector[10] means; // Means for variables
}

model {
  // Priors
  means ~ normal(0, 10);

  // Loop over constraints
  for (i in 1:N_constraints) {
    target += normal_lcdf(means[var2[i]] - means[var1[i]] | 0, sigma);
  }
}

generated quantities{
  vector[9] diff; // diffrence between means
  for (n in 1:9){
    diff[n] = means[n+1] - means[1];
  }
}

