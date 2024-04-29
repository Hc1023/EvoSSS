data {
  int<lower=0> N;           // total population
  int I0;                   // initial number of infected individuals
  real<lower=0> S0;         // initial number of susceptible individuals
  int<lower=0> T;           // number of time steps
  int<lower=0> cases[T];    // observed new cases per day
}

parameters {
  real<lower=0> beta;       // transmission rate
  real<lower=0> gamma;      // recovery rate
}

model {
  vector[T] S;
  vector[T] I;
  vector[T] R;

  // Initial conditions
  S[1] = S0;
  I[1] = I0;
  R[1] = 0;

  for (t in 1:(T-1)) {
    S[t+1] = S[t] - beta * S[t] * I[t] / N;
    I[t+1] = I[t] + beta * S[t] * I[t] / N - gamma * I[t];
    R[t+1] = R[t] + gamma * I[t];

    // likelihood of observed data
    cases[t+1] ~ poisson(beta * S[t] * I[t] / N);
  }

  // Priors
  beta ~ normal(0.4, 0.1);
  gamma ~ normal(0.1, 0.05);
}
