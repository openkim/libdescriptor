#ifndef MATHS_HPP
#define MATHS_HPP

#include <cmath>
#include "gl_quad.hpp"
#include "spherical_harmonics.hpp"
#include "radial_basis_functions.hpp"

double spherical_in(double n, double x) {
    //Modified spherical Bessel function of the first kind
    return std::cyl_bessel_i(n + 0.5, x) * sqrt(M_PI / (2 * x));
}
double spherical_in(int n, double x){
    return spherical_in(static_cast<double>(n), x);
}
#endif //MATHS_HPP
