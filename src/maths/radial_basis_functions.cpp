#include <cmath>
#include <Eigen/Core>
#include "Eigen/src/Core/Map.h"
#include "Eigen/Eigenvalues"
#include "maths/bessel_functions.hpp"
#include "maths/radial_basis_functions.hpp"

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
}

//void bessel_basis(int n_max, double rc, int r_size, double *r, double *r_basis) {
//    Eigen::VectorXd r_vec(r_size);
//    for (int i = 0; i < r_size; i++) {
//        r_vec(i) = r[i];
//    }
//
//    Eigen::MatrixXd u_all = Eigen::MatrixXd::Zero(n_max + 1, n_max + 2);
//    spherical_jn_zeros(n_max, u_all.data());
//
//    Eigen::MatrixXd coeff_a = Eigen::MatrixXd::Zero(n_max+1, n_max+1);
//    Eigen::MatrixXd coeff_b = Eigen::MatrixXd::Zero(n_max+1, n_max+1);
//
//
//    for (int l = 0; l < n_max + 1; l++) {
//        Eigen::VectorXd u = u_all.row(l);
//        for (int n = 0; n < n_max - l + 1; n++) {
//            double c = std::sqrt(2 / (rc * rc * rc * (u(n) * u(n) + u(n+1) * u(n+1))));
//            coeff_a(n, l) = u(n+1) / spherical_jn(l+1, u(n)) * c;
//            coeff_b(n, l) = -u(n) / spherical_jn(l+1, u(n+1)) * c;
//        }
//    }
//
//    Eigen::VectorXd e = Eigen::VectorXd::Zero(n_max+1);
//    Eigen::VectorXd d = Eigen::VectorXd::Ones(n_max+1);
//    Eigen::MatrixXd gnl_mat = Eigen::MatrixXd::Zero(r_size, (n_max+1) * (n_max+1));
//    Eigen::VectorXd arg = r_vec / rc;
//
//    for (int l = 0; l < n_max + 1; l++) {
//        Eigen::VectorXd u = u_all.row(l);
//        for (int n = 1; n < n_max-l+1; n++) {
//            e(n) = (u(n-1) * u(n-1) * u(n+1) * u(n+1)) / ((u(n-1) * u(n-1) + u(n) * u(n)) * (u(n) * u(n) + u(n+1) * u(n+1)));
//            d(n) = 1. - e(n) / d(n-1);
//        }
//
//        for (int n = 0; n <= n_max-l; n++) {
//            for (int a = 0; a < r_size; a++) {
//                gnl_mat(a, n*(n_max + 1) + l) = coeff_a(n, l) * spherical_jn(l, u(n) * arg(a)) + coeff_b(n, l) * spherical_jn(l, u(n+1) * arg(a));
//            }
//        }
//
//        for (int n = 1; n <= n_max-l; n++) {
//            for (int a = 0; a < r_size; a++) {
//                gnl_mat(a, n*(n_max + 1) + l) = (gnl_mat(a, n*(n_max + 1) + l) + std::sqrt(1. - d(n)) * gnl_mat(a, (n-1)*(n_max + 1) + l)) / std::sqrt(d(n));
//            }
//        }
//    }
//
//    // Copy to the output pointer
//    for (int i = 0; i < r_size; i++) {
//        for (int j = 0; j < (n_max+1) * (n_max+1); j++) {
//            r_basis[i * (n_max+1) * (n_max+1) + j] = gnl_mat(i, j);
//        }
//    }
//}
//

void bessel_basis(int n_max, double rc, int r_size, double *r, int r_basis_size, double *r_basis) {
    std::vector<double> r_vec(r, r + r_size);

    int u_all_size = (n_max + 2) * (n_max + 1);
    std::vector<double> u_all(u_all_size, 0.0);
//    spherical_jn_zeros(n_max, u_all.data());

    int coeff_size = (n_max + 1) * (n_max + 1);
    std::vector<double> coeff_a(coeff_size, 0.0);
    std::vector<double> coeff_b(coeff_size, 0.0);

    for (int l = 0; l < n_max + 1; l++) {
        for (int n = 0; n < n_max - l + 1; n++) {
            double u_n  = u_all[l + (n_max + 1) * (n)];
            double u_n1 = u_all[l + (n_max + 1) * (n  + 1)];
            double c = std::sqrt(2 / (rc * rc * rc * (u_n * u_n + u_n1 * u_n1)));
            coeff_a[n * (n_max + 1) + l] = u_n1 / spherical_jn(l + 1, u_n) * c;
            coeff_b[n * (n_max + 1) + l] = -u_n / spherical_jn(l + 1, u_n1) * c;
        }
    }

    std::vector<double> e(n_max + 1, 0.0);
    std::vector<double> d(n_max + 1, 1.0);
    std::vector<double> gnl_mat(r_size * (n_max + 1) * (n_max + 1), 0.0);

    for (int l = 0; l < n_max+1; l++) {
        for (int n = 1; n < n_max - l+1; n++) {
            double u_n   = u_all[l + (n_max + 1) * (n)];
            double u_n1  = u_all[l + (n_max + 1) * (n + 1)];
            double u_n_1 = u_all[l + (n_max + 1) * (n - 1)];
            e[n] = (u_n_1 * u_n_1 * u_n1 * u_n1) / ((u_n_1 * u_n_1 + u_n * u_n) * (u_n * u_n + u_n1 * u_n1));
            d[n] = 1.0 - e[n] / d[n - 1];
        }

        for (int n = 0; n < n_max - l + 1; n++) {
            for (int a = 0; a < r_size; a++) {
                double u_n   = u_all[l + (n_max + 1) * (n)];
                double u_n1  = u_all[l + (n_max + 1) * (n + 1)];
                gnl_mat[a * (n_max + 1) * (n_max + 1) + n * (n_max + 1) + l] = coeff_a[n * (n_max + 1) + l] * spherical_jn(l, u_n * r_vec[a] / rc) + coeff_b[n * (n_max + 1) + l] * spherical_jn(l, u_n1 * r_vec[a] / rc);
            }
        }

        for (int n = 1; n < n_max - l + 1; n++) {
            for (int a = 0; a < r_size; a++) {
                gnl_mat[a * (n_max + 1) * (n_max + 1) + n * (n_max + 1) + l] = (gnl_mat[a * (n_max + 1) * (n_max + 1) + n * (n_max + 1) + l] + std::sqrt(1.0 - d[n]) * gnl_mat[a * (n_max + 1) * (n_max + 1) + (n - 1) * (n_max + 1) + l]) / std::sqrt(d[n]);
            }
        }
    }

    for (int i = 0; i < r_basis_size; i++) {
        r_basis[i] = gnl_mat[i];
    }
}