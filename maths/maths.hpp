#ifndef MATHS_HPP
#define MATHS_HPP

#include <cmath>
#include "gl_quad.hpp"
#include "spherical_harmonics.hpp"
#include "radial_basis_functions.hpp"
#include "modified_bessel_i.hpp"

double spherical_in(double n, double x) {
    //Modified spherical Bessel function of the first kind
    return bessel_I(n + 0.5, x) * sqrt(M_PI / (2 * x));
}
double spherical_in(int n, double x){
    return spherical_in(static_cast<double>(n), x);
}

//double spherical_in_cpp(int n, double x) {
//    //Modified spherical Bessel function of the first kind
//    return std::cyl_bessel_i(static_cast<double>(n) + 0.5, x) * sqrt(M_PI / (2 * x));
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
//#include "maths.hpp"
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
#endif //MATHS_HPP
