functions {
  matrix update_fun(vector pars, matrix states_old, vector N, int poolday){
    // Function to update disease states
    // poolday : number of hostpots
    vector[poolday] S = to_vector(states_old[, 1]);  // Susceptible individuals
    vector[poolday] I1 = to_vector(states_old[, 2]); // Infected individuals, lineage 1
    vector[poolday] I2 = to_vector(states_old[, 3]); // Infected individuals, lineage 2

    real beta1 = pars[1]; // Transmission rate for lineage 1
    real beta2 = pars[2]; // Transmission rate for lineage 2
    real gamma = pars[3]; // Recovery rate

    matrix[poolday, 5] updates;

    for (i in 1:poolday){
      updates[i, 1] = S[i] - beta1 * S[i] * I1[i] / N[i] - beta2 * S[i] * I2[i] / N[i]; // Update for Susceptible
      updates[i, 2] = I1[i] + beta1 * S[i] * I1[i] / N[i] - gamma * I1[i];              // Update for Infected lineage 1
      updates[i, 3] = I2[i] + beta2 * S[i] * I2[i] / N[i] - gamma * I2[i];              // Update for Infected lineage 2
      updates[i, 4] = beta1 * S[i] * I1[i] / N[i];                                      // Onset of new cases for lineage 1
      updates[i, 5] = beta2 * S[i] * I2[i] / N[i];                                      // Onset of new cases for lineage 2
    }

    return updates;
  }

  matrix simu(matrix seed_mat_I1, matrix seed_mat_I2, vector N, int poolday, vector pars){
    // Initial conditions
    // seeds for number of poolday
    matrix[poolday, 2] Onsets_mat;
    matrix[poolday, 5] states_old;

    // Initialize states_old
    states_old[1, 1] = N[1] - seed_mat_I1[1,1] - seed_mat_I2[1,1];
    states_old[1, 2] = seed_mat_I1[1,1];
    states_old[1, 3] = seed_mat_I2[1,1];  
    
    for (i in 2:poolday){
      states_old[i, 1] = N[i];
      states_old[i, 2] = 0;
      states_old[i, 3] = 0;
    }

    // Initialize Onsets matrix
    Onsets_mat[1, 1] = 0;
    Onsets_mat[1, 2] = 0;

    // Simulate for poolday
    for (t in 2:poolday){
      // Update states
      states_old = update_fun(pars, states_old, N, poolday);

      // migrate seeds in t th day at t th hotspot
      states_old[t, 2] += seed_mat_I1[t, t];
      states_old[t, 3] += seed_mat_I2[t, t];

      // pool onsets of hotspots
      Onsets_mat[t, 1] = sum(states_old[, 4]);
      Onsets_mat[t, 2] = sum(states_old[, 5]);
    }

    return Onsets_mat;
  }
}


data {
  int<lower = 1> poolday;
  real<lower = 1> contact_init;
  // expected_matrix : observed_matrix - last onsets 
  int expected_matrix[poolday, 2]; 
  vector[3] pars;
  matrix[poolday, poolday] seed_mat_I1;
  matrix[poolday, poolday] seed_mat_I2;
  vector[poolday] seed_vec;
}

parameters {
  real<lower = 1> contact; // potential susceptible contacts
}

model {
  vector[poolday] N;
  matrix[poolday, 2] simu_mat; 

  N = seed_vec * contact + 1; // >=1
  simu_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars); 
      
  for (i in 2 : poolday) {
    // likelihood of observed data
    expected_matrix[i, 1] ~ poisson(simu_mat[i, 1]);
    expected_matrix[i, 2] ~ poisson(simu_mat[i, 2]);
  }

  // Priors
  contact ~ normal(contact_init, 5);
}


