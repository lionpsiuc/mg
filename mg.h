#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Grid structure holding data for a single level.
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
 * @brief Right-hand side function for the Poisson problem.
 *
 * @param[in] x1 First coordinate.
 * @param[in] x2 Second coordinate.
 *
 * @return double Value of the function at the given coordinates.
 */
double f(double x1, double x2);

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
void poisson(const double* x, double* y, int N, double h);

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
void rhs(double* b, int N, double h, double (*f)(double, double));

/**
 * @brief Performs the weighted Jacobi smoothing iterations.
 *
 * Applies nu iterations of the weighted Jacobi method.
 *
 * @param[in,out] level Pointer to the grid structure for the current level.
 * @param[in]     omega Weighting factor.
 * @param[in]     nu    Number of smoothing iterations.
 */
void smooth(grid* level, double omega, int nu);

/**
 * @brief Computes the residual r=b-Ax.
 *
 * @param[in]  level Pointer to the grid structure.
 * @param[out] r     Output residual vector.
 */
void residual(grid* level);

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
                 int N_fine);
