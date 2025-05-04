#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "common.h"

/**
 * @brief Creates a dynamically allocated square matrix of doubles.
 *
 * @param[in] n Size of the square matrix.
 *
 * @return A pointer to the newly allocated matrix.
 */
double** create_matrix(int n);

/**
 * @brief Frees all memory associated with a dynamically allocated matrix.
 *
 * @param[in] mat Pointer to the matrix to be freed.
 * @param[in] n   Size of the square matrix.
 */
void free_matrix(double** mat, int n);

/**
 * @brief Right-hand side function for the Poisson problem.
 *
 * @param[in] x1 First coordinate.
 * @param[in] x2 Second coordinate.
 *
 * @return Value of the function at the given coordinates.
 */
double f(double x1, double x2);

/**
 * @brief Analytical solution for the function used in the Poisson problem.
 *
 * @param[in] x1 First coordinate.
 * @param[in] x2 Second coordinate.
 *
 * @return Analytical solution at the given coordinates.
 */
double f_analytical(double x1, double x2);

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
            double omega, int nu);

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
double residual(double** x, double** b, double** r, int n, double h2);

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
void restriction(double** r_fine, double** b_coarse, int n_fine);

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
void prolongate(double** x_coarse, double** x_fine, int n_fine);

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
void gauss(double** x, double** b, int n, double h2);

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
            double omega, int nu, int l, int lmax);

#endif // MULTIGRID_H
