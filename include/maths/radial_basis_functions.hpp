#ifndef RADIAL_BASIS_FUNCTIONS_HPP
#define RADIAL_BASIS_FUNCTIONS_HPP
#include <Eigen/Core>

Eigen::MatrixXd matrixSqrt(const Eigen::MatrixXd& A);

// polynomial basis
void polynomial_basis(int n_max, double cutoff, int r_size, double *r, double *r_basis);

void bessel_basis(int n_max, double rc, int r_size, double *r, double *r_basis);

#endif // MATHS_HPP