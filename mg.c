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
