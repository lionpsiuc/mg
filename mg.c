#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mg.h"

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
 * @brief Applies the discrete Poisson operator (i.e., a matrix-vector product).
 *
 * Implements the matrix-vector product y=Ax where A is the discrete Laplacian
 * matrix using the standard five-point stencil. This function avoids explicitly
 * storing the matrix A.
 *
 * @param[in]  x Input vector.
 * @param[out] y Output vector (i.e., Ax).
 * @param[in]  N Number of grid points in each dimension.
 * @param[in]  h Grid spacing.
 */
void poisson(const double* x, double* y, int N, int h) {
  double h2_inv = 1.0 / (h * h);
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      int index = j * N + i; // Needed for row-major ordering

      // For the diagonal term, we multiply by 4/h^2
      y[index] = 4.0 * h2_inv * x[index];

      // Connect to left neighbor, if not at the left edge
      if (i > 0) {
        y[index] -= h2_inv * x[index - 1];
      }

      // Connect to right neighbor, if not at the right edge
      if (i < N - 1) {
        y[index] -= h2_inv * x[index + 1];
      }

      // Connect to top neighbor, if not at the top edge
      if (j > 0) {
        y[index] -= h2_inv * x[index - N];
      }

      // Connect to bottom neighbor, if not at bottom edge
      if (j < N - 1) {
        y[index] -= h2_inv * x[index + N];
      }
    }
  }
}

/**
 * @brief Constructs the right-hand side vector for the linear system.
 *
 * Populates the vector b with values of the function f evaluated at grid
 * points.
 *
 * @param[out] b Vector to be filled with right-hand side values.
 * @param[in]  N Number of grid points in each dimension.
 * @param[in]  h Grid spacing.
 * @param[in]  f Pointer to function evaluating the right-hand side.
 */
void rhs(double* b, int N, double h, double (*f)(double, double)) {
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      int    index = j * N + i;
      double x1    = (i + 1) * h; // The x-coordinate of the grid point
      double x2    = (j + 1) * h; // The y-coordinate of the grid point
      b[index]     = f(x1, x2);   // Evaluate our function at this point
    }
  }
}

/**
 * @brief Performs the weighted Jacobi smoothing iterations.
 *
 * Applies nu iterations of the weighted Jacobi method.
 *
 * @param[in,out] level Pointer to the grid structure for the current level.
 * @param[in]     omega Weighting factor.
 * @param[in]     nu    Number of smoothing iterations.
 */
void smooth(grid* level, double omega, int nu) {

  // Extract data from the structure for convenience
  int           N     = level->N;
  int           size  = level->size;
  double        h     = level->h;
  double        D_inv = (h * h) / 4.0;  // Inverse of diagonal
  double*       x     = level->x;       // Pointer to current solution guess
  const double* b     = level->b;       // Pointer to right-hand side
  double*       r     = level->r;       // Pointer to storage for residual
  double*       Ax    = level->Ax_temp; // Pointer to storage for Ax

  for (int k = 0; k < nu; ++k) {

    // 1. Calculate Ax and store in Ax_temp
    poisson(x, Ax, N, h);

    // 2. Calculate residual r=b-Ax and store in r
    for (int i = 0; i < size; ++i) {
      r[i] = b[i] - Ax[i];
    }

    // 3. Update solution
    for (int i = 0; i < size; ++i) {
      x[i] += omega * D_inv * r[i];
    }
  }
}

/**
 * @brief Computes the residual r=b-Ax.
 *
 * @param[in]  level Pointer to the grid structure.
 * @param[out] r     Output residual vector.
 */
void residual(grid* level) {
  int           N    = level->N;
  int           size = level->size;
  double        h    = level->h;
  const double* x    = level->x;
  const double* b    = level->b;
  double*       r    = level->r;
  double*       Ax   = level->Ax_temp;

  // Calculate Ax
  poisson(x, Ax, N, h);

  // Calculate r=b-Ax
  for (int i = 0; i < size; ++i) {
    r[i] = b[i] - Ax[i];
  }
}

/**
 * @brief Restricts the residual from a fine grid to a coarse grid using
 *        four-point averaging.
 *
 * Each coarse grid value is the average of the four corresponding fine grid
 * values.
 *
 * @param[out] r_coarse Residual vector on the coarse grid.
 * @param[in]  r_fine   Residual vector on the fine grid.
 * @param[in]  N_coarse Dimension of the coarse grid.
 * @param[in]  N_fine   Dimension of the fine grid.
 */
void restriction(double* r_coarse, const double* r_fine, int N_coarse,
                 int N_fine) {

  // Loop over coarse gird points
  for (int yc = 0; yc < N_coarse; ++jc) {   // Rows
    for (int xc = 0; xc < N_coarse; ++xc) { // Columns
      int coarse_idx = yc * N_coarse + xc;  // 1D index for coarse grid point

      // Indices of the four corresponding fine grid points
      int yf_top   = 2 * yc;     // Row index of top two fine points
      int yf_bot   = 2 * yc + 1; // Row index of bottom two fine points
      int xf_left  = 2 * xc;     // Column index of left two fine points
      int xf_right = 2 * xc + 1; // Column index of right two fine points

      // Check bounds
      if (yf_top >= N_fine || yf_bot >= N_fine || xf_left >= N_fine ||
          xf_right >= N_fine) {
        fprintf(stderr,
                "Error: Fine grid index out of bounds during restriction\n");
        r_coarse[coarse_idx] = 0.0; // Assign a default value
        continue;                   // Skip to the next coarse grid point
      }

      // Calculate indices in the 1D fine array
      int fine_idx_tl = yf_top * N_fine + xf_left;  // Top-left
      int fine_idx_tr = yf_top * N_fine + xf_right; // Top-Right
      int fine_idx_bl = yf_bot * N_fine + xf_left;  // Bottom-Left
      int fine_idx_br = yf_bot * N_fine + xf_right; // Bottom-Right

      // Average the values of the four fine grid points
      r_coarse[coarse_idx] = 0.25 * (r_fine[fine_idx_tl] + r_fine[fine_idx_tr] +
                                     r_fine[fine_idx_bl] + r_fine[fine_idx_br]);
    }
  }
}
