data {
  int<lower=0> N;           // total population
  int I10;                   // initial number of infected individuals of strain1
  int I20;                   // initial number of infected individuals of strain2
  real<lower=0> S0;         // initial number of susceptible individuals
  int<lower=0> T;           // number of time steps
  int<lower=0> cases[T];    // observed new cases per day
}

parameters {
  real<lower=0> beta1;       // transmission rate of strain1
  real<lower=0> beta2;       // transmission rate of strain2
  real<lower=0> gamma;       // recovery/removed rate
}

model {
  vector[T] S;
  vector[T] I1;
  vector[T] I2;
  vector[T] R;
  real total_cases1 = 0;
  real total_cases2 = 0;
  real expected_cases[T];

  // Initial conditions
  S[1] = S0;
  I1[1] = I10;
  I2[1] = I20;
  R[1] = 0;

  for (t in 1:(T-1)) {
    S[t+1] = S[t] - beta1 * S[t] * I1[t] / N - beta2 * S[t] * I2[t] / N;
    I1[t+1] = I1[t] + beta1 * S[t] * I1[t] / N - gamma * I1[t];
    I2[t+1] = I2[t] + beta2 * S[t] * I2[t] / N - gamma * I2[t];
    R[t+1] = R[t] + gamma * I1[t] + gamma * I2[t];


    expected_cases[t+1] = beta1 * S[t] * I1[t] / N + beta2 * S[t] * I2[t] / N;
    total_cases1 += beta1 * S[t] * I1[t] / N;
    total_cases2 += beta2 * S[t] * I2[t] / N;

    // likelihood of observed data
    cases[t+1] ~ poisson(expected_cases[t+1]);

  }

  // Constraint to maintain the total cases of strain 1 to strain 2 in the ratio 3:7
  target += normal_lpdf(total_cases1 / (total_cases1 + total_cases2) | 0.3, 0.005);

  // Priors
  beta1 ~ normal(0.4, 0.1);
  beta2 ~ normal(0.4, 0.1);
  gamma ~ normal(0.1, 0.05);
}
