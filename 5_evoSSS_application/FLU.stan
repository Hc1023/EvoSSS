functions {
  matrix simu(matrix seed_mat_I1, matrix seed_mat_I2, matrix seed_mat_I3, matrix seed_mat_I4, vector N, int poolday, real gamma, int nday, real beta1, real beta2, real beta3, real beta4){
    // Initial conditions
    // seeds for number of poolday
    matrix[nday, 4] Onsets_mat;

    vector[poolday] S;
    vector[poolday] I1;
    vector[poolday] I2;
    vector[poolday] I3;
    vector[poolday] I4;

    //  update states_new
    vector[poolday] S_new;
    vector[poolday] I1_new;
    vector[poolday] I2_new;
    vector[poolday] I3_new;
    vector[poolday] I4_new;

    // Initialize states_old
    S[1] = N[1] - seed_mat_I1[1,1] - seed_mat_I2[1,1] - seed_mat_I3[1,1] - seed_mat_I4[1,1];
    I1[1] = seed_mat_I1[1,1];
    I2[1] = seed_mat_I2[1,1];  
    I3[1] = seed_mat_I3[1,1];  
    I4[1] = seed_mat_I4[1,1];  
    
    for (i in 2:poolday){
      S[i] = N[i];
      I1[i] = 0;
      I2[i] = 0;
      I3[i] = 0;
      I4[i] = 0;
    }

    // Initialize Onsets matrix
    Onsets_mat[1, 1] = 0;
    Onsets_mat[1, 2] = 0;
    Onsets_mat[1, 3] = 0;
    Onsets_mat[1, 4] = 0;

    // Simulate for poolday
    for (t in 2:nday){
      // Update states
      S_new = S - beta1 * S .* I1 ./ N - beta2 * S .* I2 ./ N - beta3 * S .* I3 ./ N  - beta4 * S .* I4 ./ N; // Update for Susceptible 
      I1_new = I1 + beta1 * S .* I1 ./ N - gamma * I1;       // Update for Infected variant 1
      I2_new = I2 + beta2 * S .* I2 ./ N - gamma * I2;       // Update for Infected variant 2
      I3_new = I3 + beta3 * S .* I3 ./ N - gamma * I3;       // Update for Infected variant 3
      I4_new = I4 + beta4 * S .* I4 ./ N - gamma * I4;       // Update for Infected variant 4


      // pool onsets of hotspots
      Onsets_mat[t, 1] = sum(beta1 * S .* I1 ./ N); 
      Onsets_mat[t, 2] = sum(beta2 * S .* I2 ./ N); 
      Onsets_mat[t, 3] = sum(beta3 * S .* I3 ./ N); 
      Onsets_mat[t, 4] = sum(beta4 * S .* I4 ./ N); 

      // migrate seeds in t th day at t th hotspot
      if (t <= poolday) {
        I1_new[t] += seed_mat_I1[t, t];
        I2_new[t] += seed_mat_I2[t, t];
        I3_new[t] += seed_mat_I3[t, t];
        I4_new[t] += seed_mat_I4[t, t];
      }

      S = S_new;
      I1 = I1_new;
      I2 = I2_new;
      I3 = I3_new;
      I4 = I4_new;
    }

    return Onsets_mat;
  }
}

data {
  int<lower = 1> poolday;
  int<lower = 1> nday;
  // expected_matrix : observed_matrix - last onsets 
  int expected_total[4];
  matrix[poolday, poolday] seed_mat_I1;
  matrix[poolday, poolday] seed_mat_I2;
  matrix[poolday, poolday] seed_mat_I3;
  matrix[poolday, poolday] seed_mat_I4;
  vector[poolday] seed_vec;
  real<lower = 0> gamma;
  vector[5] pars_last;
}

parameters {
  real<lower = 1> contact; // potential susceptible contacts
  real<lower = 0, upper = 1> beta1;
  real<lower = 0, upper = 1> beta2;
  real<lower = 0, upper = 1> beta3;
  real<lower = 0, upper = 1> beta4;
}


model {
  vector[poolday] N;
  matrix[nday, 4] simu_mat;

  N = seed_vec * contact + 1; // >=1
  simu_mat = simu(seed_mat_I1, seed_mat_I2, seed_mat_I3, seed_mat_I4, N, poolday, gamma, nday, beta1, beta2, beta3, beta4); 
  expected_total[1] ~ poisson(sum(simu_mat[, 1]));
  expected_total[2] ~ poisson(sum(simu_mat[, 2]));
  expected_total[3] ~ poisson(sum(simu_mat[, 3]));
  expected_total[4] ~ poisson(sum(simu_mat[, 4]));

  // Priors
  contact ~ normal(pars_last[1], 10);
  beta1 ~ normal(pars_last[2], 0.01);
  beta2 ~ normal(pars_last[3], 0.01);
  beta3 ~ normal(pars_last[4], 0.01);
  beta4 ~ normal(pars_last[5], 0.01);

}


