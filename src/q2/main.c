#include "../../include/common.h"
#include "../../include/multigrid.h"

int main() {
  int    nu              = 4;
  double omega           = 2.0 / 3.0;
  double target_residual = 1e-7;
  int    max_v_cycles    = 10000;

  // Grid sizes to test
  int N_values[] = {16, 32, 64, 128, 256};
  int num_N      = sizeof(N_values) / sizeof(N_values[0]);

  // File handling
  FILE* outfile = fopen("data/comparison.csv", "w");
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
