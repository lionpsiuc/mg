#include "../../include/common.h"
#include "../../include/multigrid.h"

int main() {
  int    N_interior      = 128;
  int    n               = N_interior + 2;
  double h               = 1.0 / (N_interior + 1.0);
  double h2              = h * h;
  int    nu              = 4;
  double omega           = 2.0 / 3.0;
  double target_residual = 1e-7;
  int    max_v_cycles    = 10000;

  // File handling
  FILE* summary = fopen("data/summary.csv", "w");
  if (summary == NULL) {
    perror("Error opening summary.csv");
    return 1;
  }
  FILE* residuals = fopen("data/residuals.csv", "w");
  if (residuals == NULL) {
    perror("Error opening residuals.csv");
    fclose(summary); // Close the already opened file
    return 1;
  }
  fprintf(summary, "Levels,Cycles,Runtime,FinalResidual,CoarseSolves\n");
  fprintf(residuals, "Levels,CycleIteration,ResidualNorm\n");

  // Allocate matrices for the finest grid
  double** x    = create_matrix(n); // Solution array
  double** b    = create_matrix(n); // Right-hand side array
  double** r    = create_matrix(n); // Residual array
  double** temp = create_matrix(n); // Temporary array for smoothing

  // Check for allocation success
  if (!x || !b || !r || !temp || !x[0] || !b[0] || !r[0] || !temp[0]) {
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

  // Calculate maximum possible levels
  int max_possible_lmax = 0;
  int current_N = N_interior; // Start with the number of interior points
  while (current_N >= 2) {
    max_possible_lmax++; // Count this level
    if (current_N == 2) {
      break;
    }
    current_N /= 2;
  }
  printf("lmax for N = %d is %d\n\n", N_interior, max_possible_lmax);

  // First question
  int start_lmax = 2;
  if (max_possible_lmax <
      start_lmax) { // Test to make sure we can carry on with the iterations
    printf("Error: N = %d is too small\n", N_interior);
    free_matrix(x, n);
    free_matrix(b, n);
    free_matrix(r, n);
    free_matrix(temp, n);
    return 1;
  }
  printf("Testing lmax from %d to %d\n\n", start_lmax, max_possible_lmax);
  printf("%4s | %6s | %14s | %13s | %8s\n", "lmax", "Cycles", "Final Residual",
         "Coarse Solves", "Time (s)");
  printf("---------------------------------------------------------\n");

  // Iterate over different numbers of levels
  for (int lmax_test = start_lmax; lmax_test <= max_possible_lmax;
       ++lmax_test) {

    // Reset solution to zero initial guess for each lmax test
    for (int i = 1; i <= N_interior; i++) {
      for (int j = 1; j <= N_interior; j++) {
        x[i][j] = 0.0;
      }
    }

    // Reset the global coarse solve counter for this run
    coarse_level_solves = 0;

    int    v_count = 0; // Counter for cycles for this lmax
    double res_norm =
        residual(x, b, r, n, h2); // Calculate initial residual norm
    fprintf(residuals, "%d,%d,%.7e\n", lmax_test, v_count, res_norm);

    // Start the timer for this lmax run
    clock_t start_time = clock();

    // Algorithm loop; keep iterating until the residual is below the target or
    // max cycles are hit
    while (res_norm > target_residual && v_count < max_v_cycles) {
      vcycle(x, b, r, temp, n, h2, omega, nu, 0, lmax_test);
      res_norm = residual(x, b, r, n, h2);
      v_count++;
      fprintf(residuals, "%d,%d,%.7e\n", lmax_test, v_count, res_norm);
    }

    // Stop the timer
    clock_t end_time = clock();

    // Calculate elapsed time in seconds
    double elapsed_time = (double) (end_time - start_time) / CLOCKS_PER_SEC;

    // Print results row
    printf("%4d | %6d | %14.7e | %13d | %8.4f\n", lmax_test, v_count, res_norm,
           coarse_level_solves, elapsed_time);

    fprintf(summary, "%d,%d,%.4f,%.7e,%d\n", lmax_test, v_count, elapsed_time,
            res_norm, coarse_level_solves);
  }

  printf("---------------------------------------------------------\n");

  // Close files
  fclose(summary);
  fclose(residuals);

  // Free all allocated memory for the finest grid matrices
  free_matrix(x, n);
  free_matrix(b, n);
  free_matrix(r, n);
  free_matrix(temp, n);

  return 0;
}
