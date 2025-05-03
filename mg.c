#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/**
 * @brief Right-hand side function for the Poisson problem.
 *
 * @param[in] x1 First coordinate.
 * @param[in] x2 Second coordinate.
 *
 * @return double Value of the function at the given coordinates.
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
 * @return double Analytical solution at the given coordinates.
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
      for (int j = 1; j < n; j++) {
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
 * @return double The 2-norm of the residual vector.
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
