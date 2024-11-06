data {
  real V10;                   
  real V20;
  int T;           // number of time steps
  real alpha21;
  real mu;
  real r1;
  real r2;
  real ratios[3];    
  real ratios_sd[3];    
}

parameters {
  real<lower=100> K;
  real<lower=0, upper=0.9*alpha21> deltaa; 
}

model {
  vector[T] V1;
  vector[T] V2;
  real alpha12 = alpha21 - deltaa;
  real expected_ratio[T]; // 0,24,48,72

  // Initial conditions: 
  V1[1] = V10;
  V2[1] = V20;

  for (t in 1:(T-1)) {
    V1[t+1] = V1[t] + r1 * V1[t] * (1 - (V1[t] + alpha12 * V2[t]) / K) - mu * V1[t];
    V2[t+1] = V2[t] + r2 * V2[t] * (1 - (V2[t] + alpha21 * V1[t]) / K) + mu * V1[t];
  }


  target += normal_lpdf(V1[24+1]/(V1[24+1]+V2[24+1]) | ratios[1], ratios_sd[1]);
  target += normal_lpdf(V1[48+1]/(V1[48+1]+V2[48+1]) | ratios[2], ratios_sd[2]);
  target += normal_lpdf(V1[72+1]/(V1[72+1]+V2[72+1]) | ratios[3], ratios_sd[3]);

}
