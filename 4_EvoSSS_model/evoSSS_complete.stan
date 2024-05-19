functions
{
  matrix update_fun(vector pars, matrix states_old, vector N, int dim){
    // Function to update disease states
    vector[dim] S = to_vector(states_old[, 1]);  // Susceptible individuals
    vector[dim] I1 = to_vector(states_old[, 2]); // Infected individuals, lineage 1
    vector[dim] I2 = to_vector(states_old[, 3]); // Infected individuals, lineage 2

    real beta1 = pars[1]; // Transmission rate for lineage 1
    real beta2 = pars[2]; // Transmission rate for lineage 2
    real gamma = pars[3]; // Recovery rate

    matrix[dim, 5] updates;

    for (i in 1:dim){
      updates[i, 1] = S[i] - beta1 * S[i] * I1[i] / N[i] - beta2 * S[i] * I2[i] / N[i]; // Update for Susceptible
      updates[i, 2] = I1[i] + beta1 * S[i] * I1[i] / N[i] - gamma * I1[i];              // Update for Infected lineage 1
      updates[i, 3] = I2[i] + beta2 * S[i] * I2[i] / N[i] - gamma * I2[i];              // Update for Infected lineage 2
      updates[i, 4] = beta1 * S[i] * I1[i] / N[i];                                      // Onset of new cases for lineage 1
      updates[i, 5] = beta2 * S[i] * I2[i] / N[i];                                      // Onset of new cases for lineage 2
    }
    return updates;
  }

  matrix simu(matrix seed_mat_I1, matrix seed_mat_I2, vector N, int dim, vector pars){
    int ndays = dim * 2;

    // Initial conditions
    vector[dim] I1_old = to_vector(seed_mat_I1[1, ]);
    vector[dim] I2_old = to_vector(seed_mat_I2[1, ]);
    vector[dim] S_old = N - I1_old - I2_old;
    matrix[dim, 5] states_old;
    matrix[dim, 5] states_new;
    matrix[ndays, 3] Onsets_mat1;

    for (i in 1:dim){
      states_old[i, 1] = S_old[i];
      states_old[i, 2] = I1_old[i];
      states_old[i, 3] = I2_old[i];
    }

    // Initialize Onsets matrix

    Onsets_mat1[1, 1] = 1;
    Onsets_mat1[1, 2] = 0;
    Onsets_mat1[1, 3] = 0;

    // Loop through days to simulate
    for (t in 1:(ndays-1)){
      states_new = update_fun(pars, states_old, N, dim);

      // Update states
      states_old = states_new;
      for (i in 1 : dim){
        states_old[i, 2] += seed_mat_I1[t + 1, i];
        states_old[i, 3] += seed_mat_I2[t + 1, i];
      }

      Onsets_mat1[t + 1, 1] = t + 1; 
      Onsets_mat1[t + 1, 2] = sum(states_old[, 4]);
      Onsets_mat1[t + 1, 3] = sum(states_old[, 5]);

    }

    return Onsets_mat1;
  }
}


data
{
  int<lower = 1> dim;
  int<lower = 1> cycle;
  int<lower = 1> poolday;
  int<lower = 1> simu_days;
  int<lower = 1> observed_days;
  real contact_init;
  matrix[simu_days, 2] simu_matrix; 
  int observed_matrix[observed_days, 2]; 
  vector[3] pars;
  matrix[dim * 2, 3] Onsets_mat;
  matrix[dim, dim] Mobility_matrix;
}

parameters
{
  real<lower = 1> contact; // potential susceptible contacts
}


model
{
  vector[dim] Onset1;
  vector[dim] Onset2;
  vector[dim] seed_vec;
  vector[dim] p;
  vector[dim] N;
  matrix[2, dim] seed_matrix;
  matrix[dim * 2, dim] seed_mat_I1;
  matrix[dim * 2, dim] seed_mat_I2;
  matrix[dim * 2, 3] Onsets_mat_cycle = Onsets_mat;
  matrix[simu_days, 2] simu_matrix_cycle; 
  
  for (i in 1:simu_days){
    for (j in 1:2){
      simu_matrix_cycle[i,j] = simu_matrix[i,j];
    }
  }

  simu_matrix_cycle = simu_matrix;

  for (n in 1:cycle){
    for (i in 1 : dim){
      Onset1[i] = Onsets_mat_cycle[i + poolday, 2];
      Onset2[i] = Onsets_mat_cycle[i + poolday, 3];
    }

    seed_vec =  Mobility_matrix * (Onset1 + Onset2);
    p = Onset1 ./ (Onset1 + Onset2);

    for (i in 1:dim){
      seed_matrix[1, i] = seed_vec[i] * p[i];       // Scale seed_vec by p
      seed_matrix[2, i] = seed_vec[i] * (1 - p[i]); // Scale seed_vec by (1-p)
    }


    // Initialize the diagonal matrices
    for (i in 1 : dim) {
      for (j in 1 : dim) {
        seed_mat_I1[i, j] = i == j ? seed_matrix[1, i] : 0;  // Diagonal matrix part
        seed_mat_I1[dim + i, j] = 0;                         // Zero matrix part
        seed_mat_I2[i, j] = i == j ? seed_matrix[2, i] : 0;  // Assuming similar structure for seed_mat_I2
        seed_mat_I2[dim + i, j] = 0;
      }
    }

    N = (seed_vec + 1) * contact;
    Onsets_mat_cycle = simu(seed_mat_I1, seed_mat_I2, N, dim, pars); 

    for (i in 1 : (dim * 2)) {
      simu_matrix_cycle[i + poolday * n, 1] += Onsets_mat_cycle[i, 2];
      simu_matrix_cycle[i + poolday * n, 2] += Onsets_mat_cycle[i, 3];
    }

  }
      
  for (i in 1 : observed_days) {
    // likelihood of observed data
    observed_matrix[i, 1] ~ poisson(simu_matrix_cycle[i, 1]+1);
    observed_matrix[i, 2] ~ poisson(simu_matrix_cycle[i, 2]+1);
  }

  // Priors
  contact ~ normal(contact_init, 1);
}


