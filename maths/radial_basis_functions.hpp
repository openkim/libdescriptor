#include <cmath>
#include "Eigen/Core"
#include "Eigen/src/Core/Map.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include <iostream>
#include <ostream>

Eigen::MatrixXd matrixSqrt(const Eigen::MatrixXd& A) {
    // Ensure the matrix is square
    assert(A.rows() == A.cols());

    // Perform Schur decomposition
    Eigen::ComplexSchur<Eigen::MatrixXd> schur(A);

    // Get the unitary and upper triangular matrices
    Eigen::MatrixXd Q = schur.matrixU().real();
    Eigen::MatrixXd T = schur.matrixT().real();

    // Compute the square root of the diagonal of T
    for (int i = 0; i < T.rows(); ++i) {
        T(i, i) = std::sqrt(T(i, i));
    }

    // Construct the matrix square root
    Eigen::MatrixXd sqrtA = Q * T * Q.transpose();

    return sqrtA;
}

// polynomial basis
void polynomial_basis(int n_max, double cutoff, int r_size, double *r, double *r_basis) {
    //Dscribe polynomial basis. Ignoring the orthogonalization for now, as the paper says so!"""
    for (int i = 0; i < n_max; i++){
        int power = i + 3;
        for (int j = 0; j < r_size; j++){
            r_basis[i * r_size + j] = std::pow(cutoff - r[j],power);
        }
    }
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n_max, n_max);
    for (int i = 1; i < n_max + 1; i++){
        for (int j = 1; j < n_max + 1; j++){
            S(i - 1, j - 1) = (2 * std::pow((cutoff), (7 + i + j))) / (
                (5 + i + j) * (6 + i + j) * (7 + i + j)
            );
        }
    }
    S = S.inverse();
    S = matrixSqrt(S);

    Eigen::MatrixXd r_basis_map = Eigen::Map<Eigen::MatrixXd>(r_basis, r_size, n_max);
    Eigen::MatrixXd temp = r_basis_map * S;
    for (int i = 0; i < n_max; i++){
        for (int j = 0; j < r_size; j++){
            r_basis[i * r_size + j] = temp(j, i);
        }
    }
    // for (int i = 0; i < r_size; i++){
    //     std::cout << r_basis[i]  << " " << r_basis[r_size + i]<< std::endl;
    // }
    std::cout << "polynomial basis correct\n" << std::flush;
}
