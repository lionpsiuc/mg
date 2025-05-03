#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int coarse_level_solves = 0;

/**
 * @brief Creates a dynamically allocated square matrix of doubles.
 *
 * @param[in] n Size of the square matrix.
 *
 * @return A pointer to the newly allocated matrix.
 */
double** create_matrix(int n) {
  double** mat = (double**) malloc(n * sizeof(double*));
  for (int i = 0; i < n; i++) {
    mat[i] = (double*) calloc(n, sizeof(double));
  }
  return mat;
}

/**
 * @brief Frees all memory associated with a dynamically allocated matrix.
 *
 * @param[in] mat Pointer to the matrix to be freed.
 * @param[in] n   Size of the square matrix.
 */
void free_matrix(double** mat, int n) {
  for (int i = 0; i < n; i++) {
    free(mat[i]);
  }
  free(mat);
}

/**
 * @brief Right-hand side function for the Poisson problem.
 *
 * @param[in] x1 First coordinate.
 * @param[in] x2 Second coordinate.
 *
 * @return Value of the function at the given coordinates.
 */
double f(double x1, double x2) {
  return 2.0 * M_PI * M_PI * sin(M_PI * x1) * sin(M_PI * x2);
}

/**
 * @brief Analytical solution for the function used in the Poisson problem.
 *
 * @param[in] x1 First coordinate.
 * @param[in] x2 Second coordinate.
 *
 * @return Analytical solution at the given coordinates.
 */
double f_analytical(double x1, double x2) {
  return sin(M_PI * x1) * sin(M_PI * x2);
}

/**
 * @brief Performs the weighted Jacobi smoothing iterations.
 *
 * Applies nu iterations of the weighted Jacobi method.
 *
 * @param[in,out] x     Current solution array, which is updated in place.
 * @param[in]     b     Right-hand side of the linear system.
 * @param[in,out] temp  Temporary storage array for the previous iteration
 *                      values.
 * @param[in]     n     Grid size, including the boundary points.
 * @param[in]     h2    Square of the grid spacing.
 * @param[in]     omega Relaxation parameter.
 * @param[in]     nu    Number of smoothing iterations to perform.
 */
void smooth(double** x, double** b, double** temp, int n, double h2,
            double omega, int nu) {

  // Perform nu iterations of weighted Jacobi
  for (int iter = 0; iter < nu; iter++) {

    // Copy current solutions to temp
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        temp[i][j] = x[i][j];
      }
    }

    // Update interior points using a weighted average
    for (int i = 1; i < n - 1; i++) {
      for (int j = 1; j < n - 1; j++) {
        double new = 0.25 * (temp[i - 1][j] + temp[i + 1][j] + temp[i][j - 1] +
                             temp[i][j + 1] + h2 * b[i][j]);
        x[i][j]    = (1.0 - omega) * temp[i][j] + omega * new;
      }
    }
  }
}

/**
 * @brief Computes the residual r=b-Ax=b-(-Delta x) and its 2-norm.
 *
 * @param[in]  x  Current solution array.
 * @param[in]  b  Right-hand side of the linear system.
 * @param[out] r  Array to store the computed residual.
 * @param[in]  n  Grid size, including the boundary points.
 * @param[in]  h2 Square of the grid spacing.
 *
 * @return The 2-norm of the residual vector.
 */
double residual(double** x, double** b, double** r, int n, double h2) {
  double norm = 0.0;

  // Compute residual at the interior points
  for (int i = 1; i < n - 1; i++) {
    for (int j = 1; j < n - 1; j++) {

      // r=b-Ax=b-(-Delta x)=b+Delta x
      double laplacian = (x[i - 1][j] + x[i + 1][j] + x[i][j - 1] +
                          x[i][j + 1] - 4.0 * x[i][j]) /
                         h2;
      r[i][j] = b[i][j] + laplacian;

      norm += r[i][j] * r[i][j];
    }
  }

  // Set boundary residuals to zero
  for (int i = 0; i < n; i++) {
    r[i][0] = r[i][n - 1] = r[0][i] = r[n - 1][i] = 0.0;
  }

  return sqrt(norm); // 2-norm
}

/**
 * @brief Restricts the residual from a fine grid to a coarser grid using
 *        full-weighting.
 *
 * Implements the restriction operator that transfers the residual from level l
 * to level l+1 (coarser) using a full-weighting scheme. It uses a 3x3 stencil
 * with weights of 1/4 for the center point, 1/8 for the edge-adjacent points,
 * and 1/16 for the diagonally-adjacent points.
 *
 * @param[in]  r_fine   Residual array on the fine grid.
 * @param[in]  b_coarse Array to store the restricted residual on the coarse
 *                      grid.
 * @param[out] n_fine   Size of the fine grid, including boundary points.
 */
void restriction(double** r_fine, double** b_coarse, int n_fine) {
  int n_coarse = (n_fine - 2) / 2 + 2;

  // Full-weighting restriction
  for (int i = 1; i < n_coarse - 1; i++) {
    for (int j = 1; j < n_coarse - 1; j++) {
      int i_fine = 2 * i;
      int j_fine = 2 * j;

      // Averaging
      b_coarse[i][j] =
          0.0625 *
              (r_fine[i_fine - 1][j_fine - 1] + r_fine[i_fine - 1][j_fine + 1] +
               r_fine[i_fine + 1][j_fine - 1] +
               r_fine[i_fine + 1][j_fine + 1]) +
          0.125 * (r_fine[i_fine - 1][j_fine] + r_fine[i_fine][j_fine - 1] +
                   r_fine[i_fine + 1][j_fine] + r_fine[i_fine][j_fine + 1]) +
          0.25 * r_fine[i_fine][j_fine];
    }
  }

  // Set boundary values to zero
  for (int i = 0; i < n_coarse; i++) {
    b_coarse[i][0] = b_coarse[i][n_coarse - 1] = b_coarse[0][i] =
        b_coarse[n_coarse - 1][i]              = 0.0;
  }
}

/**
 * @brief Prolongs a correction from a coarse grid to a fine grid using bilinear
 *        interpolation.
 *
 * Implements the prolongation operator that transfers a correction from level
 * l+1 (coarse) to level l (fine) using bilinear interpolation. It works as
 * follows:
 *
 *   1. Copy coincident points to the fine grid.
 *   2. Interpolate them horizontally.
 *   3. Interpolate them vertically.
 *   4. Interpolate them diagonally to complete the center point.
 *
 * @param[in]     x_coarse Correction vector on the coarse grid.
 * @param[in,out] x_fine   Solution vector on the fine grid, updated by adding
 *                         the prolongated correction.
 * @param[in]     n_fine   Size of the fine grid, including boundary points.
 */
void prolongate(double** x_coarse, double** x_fine, int n_fine) {
  int n_coarse = (n_fine - 2) / 2 + 2;

  // 1. Copy coincident points to the fine grid
  for (int i_c = 0; i_c < n_coarse; i_c++) {
    for (int j_c = 0; j_c < n_coarse; j_c++) {
      int i_f = 2 * i_c;
      int j_f = 2 * j_c;
      if (i_f < n_fine && j_f < n_fine) {
        x_fine[i_f][j_f] += x_coarse[i_c][j_c];
      }
    }
  }

  // 2. Interpolate them horizontally
  for (int i_c = 0; i_c < n_coarse; i_c++) {
    for (int j_c = 0; j_c < n_coarse - 1; j_c++) {
      int i_f = 2 * i_c;
      int j_f = 2 * j_c + 1;
      if (i_f < n_fine && j_f < n_fine) {
        x_fine[i_f][j_f] += 0.5 * (x_coarse[i_c][j_c] + x_coarse[i_c][j_c + 1]);
      }
    }
  }

  // 3. Interpolate them vertically
  for (int i_c = 0; i_c < n_coarse - 1; i_c++) {
    for (int j_c = 0; j_c < n_coarse; j_c++) {
      int i_f = 2 * i_c + 1;
      int j_f = 2 * j_c;
      if (i_f < n_fine && j_f < n_fine) {
        x_fine[i_f][j_f] += 0.5 * (x_coarse[i_c][j_c] + x_coarse[i_c + 1][j_c]);
      }
    }
  }

  // 4. Interpolate them diagonally to complete the center point
  for (int i_c = 0; i_c < n_coarse - 1; i_c++) {
    for (int j_c = 0; j_c < n_coarse - 1; j_c++) {
      int i_f = 2 * i_c + 1;
      int j_f = 2 * j_c + 1;
      if (i_f < n_fine && j_f < n_fine) {
        x_fine[i_f][j_f] +=
            0.25 * (x_coarse[i_c][j_c] + x_coarse[i_c + 1][j_c] +
                    x_coarse[i_c][j_c + 1] + x_coarse[i_c + 1][j_c + 1]);
      }
    }
  }
}

/**
 * @brief Solves the linear system on the coarsest grid using the Gauss-Seidel
 *        method.
 *
 * This function implements the Gauss-Seidel iterative method to solve the
 * linear system Ax=b on the coarsest grid level of the multigrid hierarchy. The
 * iteration continues until convergence is reached or the maximum iteration
 * count is exceeded.
 *
 * @param[out] x  Solution array, initialized to zero and updated with the
 *                solution.
 * @param[in]  b  Right-hand side of the linear system.
 * @param[in]  n  Grid size, including boundary points.
 * @param[in]  h2 Square of the grid spacing.
 */
void gauss(double** x, double** b, int n, double h2) {
  int    max_iter = 5000;
  double tol      = 1e-10;

  // Initialize solution to zero
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      x[i][j] = 0.0;
    }
  }

  // Solve with Gauss-Seidel iterations
  for (int iter = 0; iter < max_iter; iter++) {
    double max_diff = 0.0;
    for (int i = 1; i < n - 1; i++) {
      for (int j = 1; j < n - 1; j++) {
        double old_val = x[i][j];
        x[i][j]        = 0.25 * (x[i - 1][j] + x[i + 1][j] + x[i][j - 1] +
                          x[i][j + 1] + h2 * b[i][j]);

        double diff = fabs(x[i][j] - old_val);
        if (diff > max_diff)
          max_diff = diff;
      }
    }
    if (max_diff < tol)
      break;
  }
}

/**
 * @brief Performs the algorithm.
 *
 * @param[in,out] x     Current solution vector, updated with corrections.
 * @param[in]     b     Right-hand side of the linear system.
 * @param[out]    r     Array for storing the residual.
 * @param[in,out] temp  Temporary array used in smoothing operations.
 * @param[in]     n     Grid size, including boundary points.
 * @param[in]     h2    Square of the grid spacing.
 * @param[in]     omega Relaxation parameter for weighted Jacobi.
 * @param[in]     nu    Number of smoothing operations.
 * @param[in]     l     Current grid level, where 0 is the finest.
 * @param[in]     lmax  Maximum grid level.
 */
void vcycle(double** x, double** b, double** r, double** temp, int n, double h2,
            double omega, int nu, int l, int lmax) {

  // 1. Smoothing using nu iterations of weighted Jacobi
  smooth(x, b, temp, n, h2, omega, nu);

  // 2. Compute the residuals
  residual(x, b, r, n, h2);

  int n_coarse = (n - 2) / 2 + 2; // Coarse grid size

  // If at coarsest level, solve
  if (l == lmax - 1) {
    gauss(temp, r, n, h2);
    coarse_level_solves++;

    // Add the correction to our solution
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        x[i][j] += temp[i][j];
      }
    }
  } else {

    // Allocate coarse grid arrays
    double** x_coarse    = create_matrix(n_coarse);
    double** b_coarse    = create_matrix(n_coarse);
    double** r_coarse    = create_matrix(n_coarse);
    double** temp_coarse = create_matrix(n_coarse);

    // 3. Restrict residual to coarse grid
    restriction(r, b_coarse, n);

    // 7. Recursive call
    vcycle(x_coarse, b_coarse, r_coarse, temp_coarse, n_coarse, h2 * 4.0, omega,
           nu, l + 1, lmax);

    // 8. Prolongate correction back to the fine grid
    prolongate(x_coarse, x, n);

    // Clean
    free_matrix(x_coarse, n_coarse);
    free_matrix(b_coarse, n_coarse);
    free_matrix(r_coarse, n_coarse);
    free_matrix(temp_coarse, n_coarse);
  }

  // 9. Smoothing again
  smooth(x, b, temp, n, h2, omega, nu);
}

int main() {
  int    nu              = 4;
  double omega           = 2.0 / 3.0;
  double target_residual = 1e-7;
  int    max_v_cycles    = 10000;

  // Grid sizes to test
  int N_values[] = {16, 32, 64, 128, 256};
  int num_N      = sizeof(N_values) / sizeof(N_values[0]);

  // File handling
  FILE* outfile = fopen("comparison.csv", "w");
  if (outfile == NULL) {
    perror("Error opening comparison.csv");
    return 1;
  }
  fprintf(outfile, "N,Config,lmax,Cycles,Runtime,FinalResidual\n");

  printf("Running Comparison for N = 16, 32, 64, 128, 256\n\n");
  printf("%3s | %13s | %4s | %6s | %8s | %14s\n", "N", "Configuration", "lmax",
         "Cycles", "Time (s)", "Final Residual");
  printf("---------------------------------------------------------------\n");

  // Loop over different finest grid sizes
  for (int n_idx = 0; n_idx < num_N; ++n_idx) {
    int    N_interior = N_values[n_idx];
    int    n          = N_interior + 2;
    double h          = 1.0 / (N_interior + 1.0);
    double h2         = h * h;

    // Allocate matrices for the current finest grid
    double** x    = create_matrix(n);
    double** b    = create_matrix(n);
    double** r    = create_matrix(n);
    double** temp = create_matrix(n);

    if (!x || !b || !r || !temp) {
      fprintf(stderr, "Error: Matrix allocation failed\n");
      if (x) {
        free_matrix(x, n);
      }
      if (b) {
        free_matrix(b, n);
      }
      if (r) {
        free_matrix(r, n);
      }
      if (temp) {
        free_matrix(temp, n);
      }
      fclose(outfile);
      return 1;
    }

    // Initialize right-hand side
    for (int i = 1; i <= N_interior; i++) {
      for (int j = 1; j <= N_interior; j++) {
        double x1 = (double) i * h;
        double x2 = (double) j * h;
        b[i][j]   = f(x1, x2);
      }
    }

    // Run the 2-level configuration
    int lmax_2lvl = 2;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        x[i][j] = 0.0;
      }
    }
    coarse_level_solves     = 0;
    int     v_count_2lvl    = 0;
    double  res_norm_2lvl   = residual(x, b, r, n, h2);
    clock_t start_time_2lvl = clock();
    while (res_norm_2lvl > target_residual && v_count_2lvl < max_v_cycles) {
      vcycle(x, b, r, temp, n, h2, omega, nu, 0, lmax_2lvl);
      res_norm_2lvl = residual(x, b, r, n, h2);
      v_count_2lvl++;
    }
    clock_t end_time_2lvl = clock();
    double  elapsed_time_2lvl =
        (double) (end_time_2lvl - start_time_2lvl) / CLOCKS_PER_SEC;

    // Print and write results
    printf("%3d | %13s | %4d | %6d | %8.4f | %14.7e\n", N_interior, "2-level",
           lmax_2lvl, v_count_2lvl, elapsed_time_2lvl, res_norm_2lvl);
    fprintf(outfile, "%d,2-level,%d,%d,%.4f,%.7e\n", N_interior, lmax_2lvl,
            v_count_2lvl, elapsed_time_2lvl, res_norm_2lvl);

    // Run the secondary configuration
    int lmax_maxlvl = 0;
    if (N_interior < 16) {
      printf("%3d | %13s | %4s | %6s | %8s | %14s\n", N_interior, "max-level",
             "N/A", "N/A", "N/A", "N/A"); // Indicate not applicable
      fprintf(outfile, "%d,max-level,%d,%d,%.4f,%.7e\n", N_interior, -1, -1,
              -1.0, -1.0);
    } else {
      // Add 1e-9 for float precision issues in log2 calculation
      lmax_maxlvl = (int) (log2((double) N_interior / 8.0) + 1.0 + 1e-9);

      // Reset solution
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          x[i][j] = 0.0;

      coarse_level_solves       = 0;
      int     v_count_maxlvl    = 0;
      double  res_norm_maxlvl   = residual(x, b, r, n, h2);
      clock_t start_time_maxlvl = clock();
      while (res_norm_maxlvl > target_residual &&
             v_count_maxlvl < max_v_cycles) {
        vcycle(x, b, r, temp, n, h2, omega, nu, 0, lmax_maxlvl);
        res_norm_maxlvl = residual(x, b, r, n, h2);
        v_count_maxlvl++;
      }
      clock_t end_time_maxlvl = clock();
      double  elapsed_time_maxlvl =
          (double) (end_time_maxlvl - start_time_maxlvl) / CLOCKS_PER_SEC;

      // Print and write results
      printf("%3s | %13s | %4d | %6d | %8.4f | %14.7e\n", "", "max-level",
             lmax_maxlvl, v_count_maxlvl, elapsed_time_maxlvl, res_norm_maxlvl);
      fprintf(outfile, "%d,max-level,%d,%d,%.4f,%.7e\n", N_interior,
              lmax_maxlvl, v_count_maxlvl, elapsed_time_maxlvl,
              res_norm_maxlvl);
    }

    // Free matrices for current N before next iteration
    free_matrix(x, n);
    free_matrix(b, n);
    free_matrix(r, n);
    free_matrix(temp, n);

  } // End loop over N_values

  printf("---------------------------------------------------------------\n");

  // Close file
  fclose(outfile);

  return 0;
}
