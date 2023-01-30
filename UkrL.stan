//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a matrix 'y' of 'nrows' and 'ncolumns'.
data {
  int <lower=0> nrows;
  int <lower=0> ncolumns;
  real y[nrows, ncolumns];
}

// The parameters accepted by the model. Our model
// accepts parameters 'beta0', 'beta2', 'bet0', 'bet1', and 'sigma' .
parameters {
  real <lower=0> psi;
  real beta0;
  real beta1;
  real bet0[nrows];
  real bet1[nrows];
  real sigma;
}

// The model to be estimated. We model the output
// 'y' to be beta distributed.
model {
  real eta[nrows, ncolumns];
  real mu[nrows, ncolumns];
  real A[nrows, ncolumns];
  real B[nrows, ncolumns];
  int  x[ncolumns];
  
  // priors
  psi ~ uniform(0,10000);
  sigma ~ uniform(0,100);
  beta0 ~ normal(0,sigma);
  beta1 ~ normal(0,sigma);
  
  for (i in 1:nrows) {
    bet0[i] ~ normal(0,sigma);
    bet1[i] ~ normal(0,sigma);
    for (j in 1:ncolumns) {
      x[j] = j-6;
      eta[i,j] = beta0 + bet0[i] + (beta1 + bet1[i]) * x[j];
      mu[i,j] = exp(eta[i,j]) / (1 + exp(eta[i,j]));
      A[i,j] = mu[i,j] * psi;
      B[i,j] = (1 - mu[i,j]) * psi;
      y[i,j] ~ beta(A[i,j],B[i,j]);
    }
  }
}
    
