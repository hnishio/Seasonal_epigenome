
data {
  int<lower=1> K;  // No. of categories
  int<lower=0> T;  // No. of time
  int<lower=1> N;  // No. of bins
  int<lower=1,upper=K> z[N,T]; // Category matrix (row: no. of bins, column: no. of time)
  vector<lower=0>[K] alpha;  // Prior of transition probabilities
}
parameters {
  simplex[K] theta[K];  // transition probabilities from K classes (theta[K]) to K classes (simplex[K])
  // K x K probabilities in total
}

model {
  
  for (k in 1:K){
    theta[k] ~ dirichlet(alpha); // transition probabilities from K classes
  }
  
  for (n in 1:N){
    for (t in 2:T){
      z[n,t] ~ categorical(theta[z[n,t-1]]); // present category is derived from the previous category
    }
  }
  
}

