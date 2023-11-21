#ifndef BESSEL_FUNCTIONS_HPP
#define BESSEL_FUNCTIONS_HPP
#define _USE_MATH_DEFINES

inline double chebev(const double *c, const int m, const double x);

// From numerical recipies in C++ 3rd edition
// omiting the besselk and derivative parts
double bessel_I(const double nu, const double x);

double bessel_J(const double nu, const double x);

double spherical_in(double n, double x);

double spherical_in(int n, double x);

double spherical_jn(double n, double x);

double spherical_jn(int n, double x);

double halleys_root(double l, double lwr_bnd, double upr_bnd);

double halleys_root(int l, double lwr_bnd, double upr_bnd);


void spherical_jn_zeros(int n_max, double * u_all );

#endif

//double spherical_in_cpp(int n, double x) {
//    //Modified spherical Bessel function of the first kind
//    return std::cyl_bessel_i(static_cast<double>(n) + 0.5, x) * sqrt(M_PI / (2 * x));
//}
//double spherical_jn_cpp(int n, double x) {
//    //Modified spherical Bessel function of the first kind
//    return std::cyl_bessel_j(static_cast<double>(n) + 0.5, x) * sqrt(M_PI / (2 * x));
//}
//double spherical_in_self(int n, double x) { // this relation is unstable over l = 5
//    if (n == 0) {
//        return std::sinh(x) / x;
//    } else if (n == 1) {
//        return std::cosh(x) / x - spherical_in_self(0, x) / x;
//    } else {
//        return spherical_in_self(n - 2, x) - (2 * n - 1) / x * spherical_in_self(n - 1, x);
//    }
//}

//#include <iostream>
//#include "bessel_functions.hpp"
//int main(){
//    // test above for n = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
//    // x = 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10
//    // n = 0
//    double x = 0.00001;
//    for (int i = 0; i < 10; i++){std::cout << spherical_in(i, x) - spherical_in_cpp(i, x) << std::endl;}
//    // n = 1
//    x = 0.0001;
//    for (int i = 0; i < 10; i++){std::cout << spherical_in(i, x) - spherical_in_cpp(i, x) << std::endl;}
//    // n = 2
//    x = 0.001;
//    for (int i = 0; i < 10; i++){std::cout << spherical_in(i, x) - spherical_in_cpp(i, x) << std::endl;}
//    // n = 3
//    x = 0.01;
//    for (int i = 0; i < 10; i++){std::cout << spherical_in(i, x) - spherical_in_cpp(i, x) << std::endl;}
//    // n = 4
//    x = 0.1;
//    for (int i = 0; i < 10; i++){std::cout << spherical_in(i, x) - spherical_in_cpp(i, x) << std::endl;}
//    // n = 5
//    x = 1;
//    for (int i = 0; i < 10; i++){std::cout << spherical_in(i, x) - spherical_in_cpp(i, x) << std::endl;}
//    // n = 6
//    x = 10;
//    for (int i = 0; i < 10; i++){std::cout << spherical_in(i, x) - spherical_in_cpp(i, x) << std::endl;}
//    return 0;
//}
// #include "bessel_functions.hpp"
// #include  <iostream>

// int main(){
//     int n_max = 3;
//     double u_all[(n_max + 2)*(n_max + 1)];

//     for (int i = 0; i < n_max + 2; i++){
//         for (int j = 0; j < n_max + 1; j++){
//             u_all[i * (n_max + 1) + j] = 0.0;
//         }
//     } 

//     spherical_jn_zeros(n_max, u_all);
//     for (int i = 0; i < n_max + 2; i++){
//         for (int j = 0; j < n_max + 1; j++){
//             std::cout << u_all[i * (n_max + 1) + j] << "  "; 
//         }
//         std::cout << "\n";
//     } 
//     return 0;
// }