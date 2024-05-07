data {
  int<lower=0> N;           // total population
  real<lower=0> S0;         // initial number of susceptible individuals
  int<lower=0> T;           // number of time steps
  int<lower=0> cases[T];    // observed new cases per day
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> I0;
}

parameters {
  real<lower=0, upper=T-70> t1_intro;  // Introduction time of strain 1
  real<lower=0, upper=T-70> t2_intro;  // Introduction time of strain 2

}

transformed parameters {
  int t1_intro_int = round(t1_intro);  // Convert to integer
  int t2_intro_int = round(t2_intro);  // Convert to integer
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
  I1[1] = 0;
  I2[1] = 0;
  R[1] = 0;

  for (t in 1:(T-1)) {
    if (t >= t1_intro && I1[t] == 0)
      I1[t] = I0;
    if (t >= t2_intro && I2[t] == 0)
      I2[t] = I0;

    S[t+1] = S[t] - beta * S[t] * I1[t] / N - beta * S[t] * I2[t] / N ;
    I1[t+1] = I1[t] + beta * S[t] * I1[t] / N - gamma * I1[t];
    I2[t+1] = I2[t] + beta * S[t] * I2[t] / N - gamma * I2[t];
    R[t+1] = R[t] + gamma * I1[t] + gamma * I2[t];

    
    // 1:68  (T-68):(T-1)
    if (t >= T-68){
        expected_cases[t+1] = beta * S[t] * I1[t] / N + beta * S[t] * I2[t] / N;
        total_cases1 += beta * S[t] * I1[t] / N;
        total_cases2 += beta * S[t] * I2[t] / N;

        // likelihood of observed data
        cases[t+1] ~ poisson(expected_cases[t+1]);
    }
  }

  // Constraint to maintain the total cases of strain 1 to strain 2 in the ratio 3:7
  target += normal_lpdf(total_cases1 / (total_cases1 + total_cases2) | 0.3, 0.005);

  // Priors
  t1_intro ~ uniform(5, 10);
  t2_intro ~ uniform(5, 10);
}
