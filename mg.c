#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
 */
void poisson(double* x, double* y, int N) {
  double h      = 1.0 / ((double) N + 1); // Grid spacing
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
 * @param[in]  f Pointer to function evaluating the right-hand side.
 */
void rhs(double* b, int N, double (*f)(double, double)) {
  double h = 1.0 / ((double) N + 1);
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
 * @brief Explain briefly.
 */
typedef struct {
  int     N;       // Grid dimension at this level
  double  h;       // Grid spacing (i.e., h=1/(N+1) for the unit square)
  int     size;    // Total number of points (i.e., N^2)
  double* x;       // Solution vector (or correction)
  double* b;       // Right-hand side vector (or restricted residual)
  double* r;       // Storage for the residual vector (i.e., r=b-Ax)
  double* Ax_temp; // Storage for the result of Ax (used in smoother/residual)
} grid;

/**
 * @brief Explain briefly.
 *
 * Further explanation, if required.
 *
 * @param[in,out] level   Explain briefly.
 * @param[in]     omega   Explain briefly.
 * @param[in]     nu      Explain briefly.
 */
void smooth(grid* level, double omega, int nu) {

  // Extract data from the structure for convenience
  int           N     = level->N;    // Grid dimension at this level
  int           size  = level->size; // It is assumed size is passed manually
  double        D_inv = (level->h * level->h) / 4.0; // Inverse of diagonal
  double*       x     = level->x; // Pointer to current solution guess
  const double* b     = level->b; // Pointer to right-hand side
  double*       r     = level->r; // Pointer to temporary storage for residual
  double*       Ax    = level->Ax_temp; // Pointer to temporary storage for Ax

  for (int k = 0; k < nu; ++k) {

    // 1. Calculate Ax and store in Ax_temp
    poisson(x, Ax, N);

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
