
data {
  int Ndate_full;
  int Ndate;
  int Nrep[Ndate];
  int sumNrep;
  vector[sumNrep] Y;
  int obsidx[Ndate];
}

parameters {
  vector<lower=0>[Ndate_full] mu;
  real<lower=0> s_mu;
  real<lower=0> s_Y;
}

model {
  // State equation of mu
  for (t in 3:Ndate_full) {
    target += normal_lpdf(mu[t]|2*mu[t-1] - mu[t-2], s_mu);
  }
  // Observation equation of Y
  for (n in 1:Nrep[1]) {                   // Nrep=c(4,4,4,...)
    target += normal_lpdf(Y[n]|mu[obsidx[1]], s_Y);
  }
  for (t in 2:Ndate) {
    for (n in 1:Nrep[t]) {  //
      target +=  normal_lpdf(Y[sum(Nrep[1:(t-1)])+n]|mu[obsidx[t]], s_Y);
    }
  }
}
