functions {
  matrix simu(matrix seed_mat_I1, matrix seed_mat_I2, vector N, int poolday, vector pars, int n){
    // Initial conditions
    // seeds for number of poolday
    matrix[n, 2] Onsets_mat;
    
    real beta1 = pars[1]; // Transmission rate for lineage 1
    real beta2 = pars[2]; // Transmission rate for lineage 2
    real gamma = pars[3]; // Recovery rate

    vector[poolday] S;
    vector[poolday] I1;
    vector[poolday] I2;

    //  update states_new
    vector[poolday] S_new;
    vector[poolday] I1_new;
    vector[poolday] I2_new;

    // Initialize states_old
    S[1] = N[1] - seed_mat_I1[1,1] - seed_mat_I2[1,1];
    I1[1] = seed_mat_I1[1,1];
    I2[1] = seed_mat_I2[1,1];  
    
    for (i in 2:poolday){
      S[i] = N[i];
      I1[i] = 0;
      I2[i] = 0;
    }

    // Initialize Onsets matrix
    Onsets_mat[1, 1] = 0;
    Onsets_mat[1, 2] = 0;

    // Simulate for poolday
    for (t in 2:n){
      // Update states
      S_new = S - beta1 * S .* I1 ./ N - beta2 * S .* I2 ./ N; // Update for Susceptible 
      I1_new = I1 + beta1 * S .* I1 ./ N - gamma * I1;       // Update for Infected lineage 1
      I2_new = I2 + beta2 * S .* I2 ./ N - gamma * I2;       // Update for Infected lineage 2

      // pool onsets of hotspots
      Onsets_mat[t, 1] = sum(beta1 * S .* I1 ./ N); 
      Onsets_mat[t, 2] = sum(beta2 * S .* I2 ./ N); 

      // migrate seeds in t th day at t th hotspot
      if (t <= poolday) {
        I1_new[t] += seed_mat_I1[t, t];
        I2_new[t] += seed_mat_I2[t, t];
      }

      S = S_new;
      I1 = I1_new;
      I2 = I2_new;
    }

    return Onsets_mat;
  }
}

data {
  int<lower = 1> poolday;
  int<lower = 1> nday;
  // expected_matrix : observed_matrix - last onsets 
  int expected_matrix[nday, 2]; 
  vector[3] pars;
  vector[poolday] Onset1;
  vector[poolday] Onset2;
  vector[2] pars_last;
}

parameters {
  real<lower = 1> contact; // potential susceptible contacts
  real<lower = 0, upper = 0.1> mobility;
}


model {
  vector[poolday] N;
  matrix[nday, 2] simu_mat;
  matrix[poolday, poolday] Mobility_matrix;
  matrix[poolday, poolday] seed_mat_I1;
  matrix[poolday, poolday] seed_mat_I2;
  vector[poolday] seed_vec;

  Mobility_matrix = diag_matrix(rep_vector(mobility, poolday));

  seed_mat_I1 = diag_matrix(Mobility_matrix * Onset1);
  seed_mat_I2 = diag_matrix(Mobility_matrix * Onset2);
  seed_vec = Mobility_matrix * (Onset1 + Onset2);

  N = seed_vec * contact + 1; // >=1
  simu_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars, nday); 
  
  for (i in 2 : nday) {
    // likelihood of observed data
    expected_matrix[i, 1] ~ poisson(simu_mat[i, 1]);
    expected_matrix[i, 2] ~ poisson(simu_mat[i, 2]);
  }

  // Priors
  mobility ~ normal(pars_last[1], 0.01);
  contact ~ normal(pars_last[2], 10);

}

