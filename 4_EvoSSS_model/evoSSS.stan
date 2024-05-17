data
{
  int<lower = 1> dim;
  int<lower = 1> cycle;
  vector[3] pars;
  matrix[dim * 2, 3] Onsets_mat;
  matrix[dim * 2, dim] seed_mat_I1;
  matrix[dim * 2, dim] seed_mat_I2;
  vector[dim] N_vec;
  matrix[dim, dim] Mobility_matrix;
}

parameters
{
  vector<lower = 1>[cycle] contact_vec; // potential susceptible contacts
}

functions
{
  // Function to update disease statesvector update_fun(vector pars, vector states_old, vector N) {
  vector[dim] S = to_vector(states_old[, 1]);  // Susceptible individuals
  vector[dim] I1 = to_vector(states_old[, 2]); // Infected individuals, lineage 1
  vector[dim] I2 = to_vector(states_old[, 3]); // Infected individuals, lineage 2

  real beta1 = pars[1]; // Transmission rate for lineage 1
  real beta2 = pars[2]; // Transmission rate for lineage 2
  real gamma = pars[3]; // Recovery rate

  matrix[dim, 5] updates; 
  for (i in 1 : dim)
  {
    updates[i, 1] = S[i] - beta1 * S[i] * I1[i] / N[i] - beta2 * S[i] * I2[i] / N[i]; // Update for Susceptible
    updates[i, 2] = I1[i] + beta1 * S[i] * I1[i] / N[i] - gamma * I1[i];              // Update for Infected lineage 1
    updates[i, 3] = I2[i] + beta2 * S[i] * I2[i] / N[i] - gamma * I2[i];              // Update for Infected lineage 2
    updates[i, 4] = beta1 * S[i] * I1[i] / N[i];                                      // Onset of new cases for lineage 1
    updates[i, 5] = beta2 * S[i] * I2[i] / N[i];                                      // Onset of new cases for lineage 2
  }
  return updates; 
}

matrix simu(matrix seed_mat_I1, matrix seed_mat_I2, vector N, int dim, vector pars)
{
  int ndays = dim * 2;

  // Initialize Onsets matrix
  matrix[ndays, 3] Onsets_mat;
  Onsets_mat[1, 1] = 1;
  Onsets_mat[1, 2] = 0;
  Onsets_mat[1, 3] = 0;

  // Initial conditions
  vector[dim] I1_old = to_vector(seed_mat_I1[1, ]);
  vector[dim] I2_old = to_vector(seed_mat_I2[1, ]);
  vector[dim] S_old = N - I1_old - I2_old;
  matrix[dim, 5] states_old;
  matrix[dim, 5] states_new;
  for (int i = 1; i <= dim; i++)
  {
    states_old[i, 1] = S_old[i];
    states_old[i, 2] = I1_old[i];
    states_old[i, 3] = I2_old[i];
  }

  // Loop through days to simulate
  for (int t = 1; t < ndays; t++)
  {
    states_new = update_fun(pars, states_old, N);

    // Update states
    states_old = states_new;
    for (int i = 1; i <= dim; i++)
    {
      states_old[i, 2] += seed_mat_I1[t + 1, i];
      states_old[i, 3] += seed_mat_I2[t + 1, i];
    }

    Onsets_mat[t + 1, 1] = t + 1
    Onsets_mat[t + 1, 2] = sum(states_old[,4])
    Onsets_mat[t + 1, 3] = sum(states_old[,5])
  }

  return Onsets_mat;
}

model
{
  int poolday = 30;
  vector[dim] Onset1;
  vector[dim] Onset2;
  matrix[dim * 3, 3] Onsets_mat_old;
  vector[dim] seed_vec;
  vector[dim] p;
  matrix[2, dim] seed_matrix;

  for (n in 1 : cycle)
  {

    Onsets_mat_old = Onsets_mat;

    for (i in(1 + poolday) : (dim + poolday))
    {
      Onset1[i] = Onsets_mat_old[i, 2];
      Onset2[i] = Onsets_mat_old[i, 3];
    }

    seed_vec = (Onset1 ' + Onset2') *Mobility_matrix
    p = Onset1./ (Onset1 + Onset2)

    for (i in 1 : dim)
    {
      seed_matrix[1, i] = seed_vec[i] * p[i];       // Scale seed_vec by p
      seed_matrix[2, i] = seed_vec[i] * (1 - p[i]); // Scale seed_vec by (1-p)
    }

    seed_mat_list = list()
    seed_mat_list[[1]] = rbind(diag(seed_matrix[1, ]), matrix(0, dim, dim))
    seed_mat_list[[2]] = rbind(diag(seed_matrix[2, ]), matrix(0, dim, dim))
    seed_mat_list[[3]] = (seed_vec + 1) * contact_vec[n] Onsets_mat = simu(seed_mat_list, dim)
    Onsets_mat[, 1] = Onsets_mat[, 1] + sum(poolday * n)
                                                                                 Onsets_mat_list [[n + 1]] = Onsets_mat
  }
}

generated quantities
{
  real onset1[T];
  real onset2[T];
  for (t in 1 : T)
  {
    onset1[t] = beta1 * S[t] * I1[t] / N;
    onset2[t] = beta2 * S[t] * I2[t] / N;
  }
}
