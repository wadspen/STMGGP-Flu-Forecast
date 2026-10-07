functions {
  
  matrix kronecker_prod(matrix A, matrix B) {
    int r_A = rows(A);
    int c_A = cols(A);
    int r_B = rows(B);
    int c_B = cols(B);
    
    matrix[r_A * r_B, c_A * c_B] Z;
    
    for (i in 1:r_A) {
      for (j in 1:c_A) {
        int row_start = (i - 1) * r_B + 1;
        int col_start = (j - 1) * c_B + 1;
        
        // Block-wise assignment
        Z[row_start:(row_start + r_B - 1), col_start:(col_start + c_B - 1)] = A[i, j] * B;
      }
    }
    return Z;
  }
  
  vector kron_vectors(vector x, vector y) {
    int size_x = num_elements(x);
    int size_y = num_elements(y);
    vector[size_x * size_y] result;
    
    for (i in 1:size_x) {
      for (j in 1:size_y) {
        result[(i - 1) * size_y + j] = x[i] * y[j];
      }
    }
    return result;
  }
  

  

  
  matrix exp_kernel(matrix dist_mat,
                    real length_scale, 
                    real outputvariance,
                    real delta) {
                                 
                                 
    int N = rows(dist_mat);
    matrix[N, N] K;
    for (i in 1:(N - 1)) {
      K[i, i] = outputvariance + delta;
      for (j in (i + 1):N) {
        K[i, j] = outputvariance * 
                    exp(-dist_mat[i, j] / 
                    length_scale);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = outputvariance + delta;
    return K;
  }
  
    matrix exp2_kernel(matrix dist_mat,
                    real length_scale, 
                    real outputvariance,
                    real delta) {
                                 
                                 
    int N = rows(dist_mat);
    matrix[N, N] K;
    for (i in 1:(N - 1)) {
      K[i, i] = outputvariance + delta;
      for (j in (i + 1):N) {
        K[i, j] = outputvariance * 
                    exp(-dist_mat[i, j]^2 / 
                    length_scale);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = outputvariance + delta;
    return K;
  }
  
  matrix mat52_kernel(matrix dist_mat,
                    real length_scale, 
                    real outputvariance,
                    real delta) {
                                 
                                 
    int N = rows(dist_mat);
    matrix[N, N] K;
    for (i in 1:(N - 1)) {
      K[i, i] = outputvariance + delta;
      for (j in (i + 1):N) {
        K[i, j] = outputvariance *
                    (1 + sqrt(5) * dist_mat[i, j] / length_scale + 
                     5*(dist_mat[i, j]^2) / (3 * length_scale^2)) * 
                    exp(-sqrt(5)*dist_mat[i, j] / 
                    length_scale);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = outputvariance + delta;
    return K;
  }
  
  matrix mat32_kernel(matrix dist_mat,
                    real length_scale, 
                    real outputvariance,
                    real delta) {
                                 
                                 
    int N = rows(dist_mat);
    matrix[N, N] K;
    for (i in 1:(N - 1)) {
      K[i, i] = outputvariance + delta;
      for (j in (i + 1):N) {
        K[i, j] = outputvariance *
                    (1 + sqrt(3) * dist_mat[i, j] / length_scale) * 
                    exp(-sqrt(3)*dist_mat[i, j] / 
                    length_scale);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = outputvariance + delta;
    return K;
  }
  
  matrix bern_kernel(matrix dist_mat, 
                     int P, 
                     real adist,
                     real delta) {
                       
    int N = rows(dist_mat);
    matrix[N, N] K;
    for (i in 1:N) {
      for (j in i:N) {
        K[i, j] = pow(1/(1 + adist^2 * dist_mat[i,j]^2), P * .5);
        K[j, i] = K[i, j];
      }
    }
    
    for (i in 1:N) {
      K[i, i] = K[i, i] + delta;
    }
    return K;
  }
  

}
data {
  int N;                                //total number of observations
  int<lower=1> nweeks;                  //number of weeks in the season
  int<lower=1> nseas;                   //number of seasons
  int<lower=1> nlocs;                   //number of locations
  vector[N] weeks;                      //season weeks
  vector[nseas] seas;                   //seasons
  vector[nlocs] locs;                   //locations
  matrix[nweeks, nweeks] tdist;         //week distance matrix
  matrix[nseas, nseas] sdist;           //season distance matrix
  matrix[nlocs, nlocs] ldist;           //population center distance matrix
  array[N] int rate;
  vector[N] pop;
  // vector[N] rate;                       //observed data
  array[nseas*nlocs] int first_ind;              //indices for data
  array[nseas*nlocs] int last_ind;               //last indices for data
  // array[nseas*nlocs] int first_rate_ind;              //indices for data
  // array[nseas*nlocs] int last_rate_ind;               //last indices for data
  real<lower=0> sig_alpha;    
}
transformed data {
  real jitter = 1e-8;
}

parameters {

  // Covariance function parameters
  real<lower=0> sigma;
  real rho;
  
  // Noise variance
  // vector<lower=0>[nseas*nlocs] sig;
  real<lower=0> sig;
  vector[nweeks*nseas*nlocs] epsilon;
  // matrix<lower=0>[nseas, nlocs] sig;
  // real<lower=0> sig;

  // Intercepts
  vector[nweeks*nseas*nlocs] eta;
  // vector[nseas] teta;
  
  real<lower=0> alpha;
  real<lower=0> gamma;
  // real<lower=0,upper=1> beta;
  // real<lower=0> delta;
  
  cholesky_factor_corr[nseas] L_A;
  
  // array[nlocs] matrix[nweeks, nseas] Z;

}
transformed parameters {

  matrix[nseas, nseas] K_A = bern_kernel(sdist, nseas - 1,
                                         alpha, jitter);
  // matrix[nseas, nseas] L_A = cholesky_decompose(K_A);
  
  matrix[nlocs, nlocs] K_C = bern_kernel(ldist, nseas - 1,
                                         gamma, 1e-6);
  matrix[nlocs, nlocs] L_C = cholesky_decompose(K_C);
                                         
  matrix[nweeks, nweeks] K_B = mat52_kernel(tdist, exp(rho), sigma, jitter);
  matrix[nweeks, nweeks] L_B = cholesky_decompose(K_B);
  
  
  // matrix[nweeks*nseas*nlocs, nweeks*nseas*nlocs] L_K = kronecker_prod(L_C, 
  //                                                     kronecker_prod(L_A, L_B));
  
  array[nlocs] matrix[nweeks, nseas] Z;

  for (c in 1:nlocs) {
    int start = (c - 1) * nweeks * nseas + 1;
    int stop  = c * nweeks * nseas;
  
    Z[c] = to_matrix(
              segment(eta, start, nweeks * nseas),
              nweeks,
              nseas
           );
  }
  
  array[nlocs] matrix[nweeks, nseas] X;
  
  for (c in 1:nlocs) {
    X[c] = L_B * Z[c] * L_A';
  }
  
  array[nlocs] matrix[nweeks, nseas] Y;

  for (w in 1:nweeks) {
    for (s in 1:nseas) {
  
      vector[nlocs] temp;
  
      for (c in 1:nlocs)
        temp[c] = X[c][w,s];
  
      temp = L_C * temp;
  
      for (c in 1:nlocs)
        Y[c][w,s] = temp[c];
    }
  }
  
  vector[nweeks*nseas*nlocs] f;

  for (c in 1:nlocs) {
    for (s in 1:nseas) {
      for (w in 1:nweeks) {
  
        int ind = w +
                  (s - 1) * nweeks +
                  (c - 1) * nweeks * nseas;
  
        f[ind] = Y[c][w,s];
      }
    }
  }

  // vector[nweeks*nseas*nlocs] f = L_K * eta;
  
  // vector[nseas] seasf = L_A * teta;
  
  // matrix[nseas, nseas] K_A = L_A * L_A';
  
  
}
model {
  // Covariance function parameters
  sigma ~ std_normal();
  rho ~ std_normal();
  // delta ~ std_normal();
  epsilon ~ normal(0, sig);


  sig ~ std_normal();
  eta ~ std_normal();
  // X ~ std_normal();
  // teta ~ std_normal();
  
  L_A ~ lkj_corr_cholesky(2.0);
  
  // beta ~ student_t(3, 0, 1);
  // alpha ~ std_normal();
  alpha ~ normal(0, sig_alpha);
  // gamma ~ student_t(3, 0, 1);
  gamma ~ cauchy(0, 1);
  

  // Data
  for (c in 1:nlocs) {
    for (s in 1:nseas) {
        
        int ind = nseas*(c - 1) + s;
        for (i in first_ind[ind]:last_ind[ind]) {
           rate[i] ~ poisson(exp(f[i] + epsilon[i]) * pop[i]);
        }
        
    }
  }
  // for (c in 1:nlocs) {
  //   for (s in 1:nseas) {
  //       
  //       int ind = nseas*(c - 1) + s;
  //       rate[first_ind[ind]:last_ind[ind]] ~ 
  //                                normal(f[first_ind[ind]:last_ind[ind]],
  //                                       sig[ind]);
  //       
  //       
  //   }
  // }
}


generated quantities {
  // vector[nweeks*nseas*nlocs] pred_f;
  array[nweeks*nseas*nlocs] int pred_f;
  for (c in 1:nlocs) {
    for (s in 1:nseas) {
      int ind = nseas*(c - 1) + s;
      int b_seas = (ind - 1)*nweeks + 1;
      int e_seas = ind*nweeks;
      // matrix[nweeks, nweeks] Sigma = diag_matrix(rep_vector(sig[ind], nweeks));
      // pred_f[b_seas:e_seas] = multi_normal_rng(f[b_seas:e_seas], Sigma);
      for (i in b_seas:e_seas) {
        pred_f[i] = poisson_rng(exp(f[i] + epsilon[i]) * pop[i]); 
      }
    }
  }
}





