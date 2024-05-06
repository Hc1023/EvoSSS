data {
  int<lower=0> N;           // total population
  int I0;                   // initial number of infected individuals of strain1
  real<lower=0> S0;         // initial number of susceptible individuals
  int<lower=0> T;           // number of time steps
  int<lower=0> cases[T];    // observed new cases per day
}

parameters {
  real<lower=0> beta;
  real<lower=0> gamma;       // recovery/removed rate
}

model {
  vector[T] S;
  vector[T] I;
  real expected_cases[T];

  // Initial conditions
  S[1] = S0;
  I[1] = I0;

  for (t in 1:(T-1)) {
    S[t+1] = S[t] - beta * S[t] * I[t] / N;
    I[t+1] = I[t] + beta * S[t] * I[t] / N - gamma * I[t];

    expected_cases[t+1] = beta * S[t] * I[t] / N;

    // likelihood of observed data
    cases[t+1] ~ poisson(expected_cases[t+1]);

  }

  // Priors
  beta ~ normal(0.4, 0.1);
  gamma ~ normal(0.1, 0.05);
}
