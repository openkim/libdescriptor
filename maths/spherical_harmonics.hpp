#ifndef SPHERICAL_HARMONICS_HPP
#define SPHERICAL_HARMONICS_HPP
//From Spherical Harmonic Lighting: The Gritty Details
// No longer the case: -and google's spherical harmonic library: https://github.com/google/spherical-harmonics-
// Based on http://phys.uri.edu/nigh/NumRec/bookfpdf/f6-8.pdf Numerical Recipes in F77

//Convention, largely follows the chemist's notation, simply putting larger numbers first
// first l, then m, then theta, then phi
// theta is the polar angle, phi is the azimuthal angle
// theta is in [0, pi], phi is in [0, 2pi], m is in [-l, l]
// Condon-Shortley phase: (-1)^m

#include <iostream>
#define _USE_MATH_DEFINES

#include <cmath>
#include <stdexcept>
#include <array>
#include <vector>

#define HARD_CODED_SH_COEFF 4
#define FACT_CACHE_SIZE 16


// factorials and double facotrials ************************************************

// The vast majority of SH evaluations will hit these precomputed values.
double factorial(int x) {
    static const double factorial_cache[FACT_CACHE_SIZE] = {1, 1, 2, 6, 24, 120, 720, 5040,
                                                            40320, 362880, 3628800, 39916800,
                                                            479001600, 6227020800,
                                                            87178291200, 1307674368000};

    if (x < FACT_CACHE_SIZE) {
        return factorial_cache[x];
    } else {
        double s = factorial_cache[FACT_CACHE_SIZE - 1];
        for (int n = FACT_CACHE_SIZE; n <= x; n++) {
            s *= n;
        }
        return s;
    }
}
double factorial(double x) {
    return factorial(static_cast<int>(x));
}

double double_factorial(int x) {
    static const double dbl_factorial_cache[FACT_CACHE_SIZE] = {1, 1, 2, 3, 8, 15, 48, 105,
                                                                384, 945, 3840, 10395, 46080,
                                                                135135, 645120, 2027025};

    if (x < FACT_CACHE_SIZE) {
        return dbl_factorial_cache[x];
    } else {
        double s = dbl_factorial_cache[FACT_CACHE_SIZE - (x % 2 == 0 ? 2 : 1)];
        double n = x;
        while (n >= FACT_CACHE_SIZE) {
            s *= n;
            n -= 2.0;
        }
        return s;
    }
}
// End factorials and double factorials ********************************************

// Misc math functions *************************************************************
// clamp the first argument to be greater than or equal to the second
// and less than or equal to the third.
double clamp(double val, double min, double max) {
    if (val < min) { val = min; }
    if (val > max) { val = max; }
    return val;
}

// Return true if the first value is within epsilon of the second value.
bool approx_equal(double actual, double expected) {
    double diff = std::abs(actual - expected);
    // 5 bits of error in mantissa (source of '32 *')
    return diff < 32 * std::numeric_limits<double>::epsilon();
}

// Return floating mod x % m.
double fmod(double x, double m) {
    return x - (m * floor(x / m));
}

std::array<double,3> sph2cart(double phi, double theta) {
    double r = sin(theta);
    return {r * cos(phi), r * sin(phi), cos(theta)};
}

std::array<double,3> cart2sph(const std::array<double,3> &r) {
    double norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    std::array<double, 3> r_norm = {r[0] / norm, r[1] / norm, r[2] / norm};
    // Explicitly clamp the z coordinate so that numeric errors don't cause it
    // to fall just outside of acos' domain.
    double const theta = acos(clamp(r_norm[2], -1.0, 1.0));
    // We don't need to divide r.y() or r.x() by sin(theta) since they are
    // both scaled by it and atan2 will handle it appropriately.
    double const phi = atan2(r_norm[1], r_norm[0]);
    return {phi, theta, norm};
}


// *********************************************************************************

// Condon-Shortley phase: (-1)^m
inline double condon_shortley(int m) {
    return (m % 2 == 0) ? 1.0 : -1.0;
}

// Evaluate Associated Legendre Polynomial P(l,m) at x
// Evaluate the associated Legendre polynomial of degree @l and order @m at
// coordinate @x. The inputs must satisfy:
// 1. l >= 0
// 2. 0 <= m <= l
// 3. -1 <= x <= 1
// See http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
//
double Plm(int l, int m, double theta) {
    // Compute Pmm(x) = (-1)^m(2m - 1)!!(1 - x^2)^(m/2), where !! is the double
    // factorial.
    double const x = cos(theta);
    double pmm = 1.0;
    // P00 is defined as 1.0, do don't evaluate Pmm unless we know m > 0
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= (-fact) * somx2;
            fact += 2.0;
        }
    }

    if (l == m) {
        // Pml is the same as Pmm so there's no lifting to higher bands needed
        return pmm;
    }

    // Compute Pmm+1(x) = x(2m + 1)Pmm(x)
    double pmm1 = x * (2 * m + 1) * pmm;
    if (l == m + 1) {
        // Pml is the same as Pmm+1 so we are done as well
        return pmm1;
    }

    // Use the last two computed bands to lift up to the next band until l is
    // reached, using the recurrence relationship:
    // Pml(x) = (x(2l - 1)Pml-1 - (l + m - 1)Pml-2) / (l - m)
    for (int n = m + 2; n <= l; n++) {
        double pmn = (x * (2 * n - 1) * pmm1 - (n + m - 1) * pmm) / (n - m);
        pmm = pmm1;
        pmm1 = pmn;
    }
    // Pmm1 at the end of the above loop is equal to Pml
    return pmm1;
}


double kron(int i, int j) {
    if (i == j) {
        return 1.0;
    } else {
        return 0.0;
    }
}


// Normalization Factor for SH Function
double K(int l, int m) {
    double temp = ((2.0 * l + 1.0) * factorial(l - m)) / (4.0 * M_PI * factorial(l + m));
    return sqrt(temp);
}

// Evaluate SH Function at (theta,phi) **********************************************

std::array<double,2> Ylmi(int l, int m, double phi, double theta) {
    if (l < 0) {
        throw std::runtime_error("l must be >= 0");
    }
    if (!(-l <= m && m <= l)) {
        throw std::runtime_error("m must be between -l and l.");
    }

    if (m < 0) {
        auto sh = Ylmi(l, -m, phi, theta);
        double phase = std::pow(-1.0, m);
        return {phase * sh[0], -phase * sh[1]};
    }
    double kml = sqrt((2.0 * l + 1) * factorial(l - m) /
                      (4.0 * M_PI * factorial(l + m)));
    double pml = Plm(l, m, theta);
    return {kml * cos(m * phi) * pml, kml * sin(m * phi) * pml};
}

double Ylm(int l, int m, double phi, double theta) {
//    double sign = (m % 2 == 0) ? 1.0 : -1.0;
    double sign = condon_shortley(m);
    auto sh = Ylmi(l, abs(m), phi, theta);
    if (m == 0) {
        return sh[0];
    }
    else if (m < 0) {
        return sign * M_SQRT2 * sh[1];
    }
    else {
        return sign * M_SQRT2 * sh[0];
    }
}

std::vector<double> Ylm_all_m(int l, double phi, double theta) {
    std::vector<double> sh(2 * l + 1);
    for (int m = -l; m <= l; m++) {
        sh[m + l] = (Ylm(l, m, phi, theta));
    }
    return sh;
}

std::vector<double> Ylmi_all_m(int l, double phi, double theta) {
    // first 2l+1 are real, second 2l+1 are imaginary
    int total_elem = (2 * l + 1);
    std::vector<double> sh(2 * total_elem);
    for (int m = -l; m <= l; m++) {
        auto sh_i = (Ylmi(l, m, phi, theta));
        sh[(m + l)] = sh_i[0];
        sh[(m + l) + total_elem] = sh_i[1];
    }
    return sh;
}

std::vector<double> Ylm_all_m_from_r(int l, const std::array<double,3> &r) {
    auto const sph = cart2sph(r);
    std::vector<double> sh(2 * l + 1);
    for (int m = -l; m <= l; m++) {
        sh[m + l] = Ylm(l, m, sph[0], sph[1]);
    }
    return sh;
}

std::vector<double> Ylmi_all_m_from_r(int l, const std::array<double,3> &r) {
    auto const sph = cart2sph(r);
    int total_elem = (2 * l + 1);
    std::vector<double> sh(2 * total_elem);
    for (int m = -l; m <= l; m++) {
        auto sh_i = (Ylmi(l, m, sph[0], sph[1]));
        sh[(m + l)] = sh_i[0];
        sh[(m + l) + total_elem] = sh_i[1];
    }
    return sh;
}

// pointer based version of above two functions
void Ylm_all_m(int l, double phi, double theta, double * sh) {
    for (int m = -l; m <= l; m++) {
        sh[m +l] = (Ylm(l, m, phi, theta));
    }
}

void Ylmi_all_m(int l, double phi, double theta, double * sh) {
    // first 2l+1 are real, second 2l+1 are imaginary
    int total_elem = (2 * l + 1);
    for (int m = -l; m <= l; m++) {
        auto sh_i = (Ylmi(l, m, phi, theta));
        sh[(m + l)] = sh_i[0];
        sh[(m + l) + total_elem] = sh_i[1];
    }
}

void Ylm_all_m_from_r(int l,  double * r, double * sh) {
    auto const sph = cart2sph({r[0], r[1], r[2]});
    for (int m = -l; m <= l; m++) {
        sh[m + l] = Ylm(l, m, sph[0], sph[1]);
    }
}

void Ylmi_all_m_from_r(int l, double * r, double * sh){
    auto const sph = cart2sph({r[0], r[1], r[2]});
    int total_elem = (2 * l + 1);
    for (int m = -l; m <= l; m++) {
        auto sh_i = (Ylmi(l, m, sph[0], sph[1]));
        sh[(m + l)] = sh_i[0];
        sh[(m + l) + total_elem] = sh_i[1];
    }
}

void Ylmi_all_l_from_r(int l_max, double * r, double * sh){
    auto const sph = cart2sph({r[0], r[1], r[2]});
    int total_elem = (l_max + 1) * (l_max + 1);
    for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
            auto sh_i = (Ylmi(l, m, sph[0], sph[1]));
            sh[(l * l + l + m)] = sh_i[0];
            sh[(l * l + l + m) + total_elem] = sh_i[1];
        }
    }
}

void Ylmi_all_l_from_r(int l_max, double * r, double * sh_real, double * sh_imag){
    auto const sph = cart2sph({r[0], r[1], r[2]});
    int total_elem = (l_max + 1) * (l_max + 1);
    for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
            auto sh_i = (Ylmi(l, m, sph[1], sph[0]));
            sh_real[(l * l + l + m)] = sh_i[0];
            sh_imag[(l * l + l + m)] = sh_i[1];
        }
    }
}

void Ylm_all_l_from_r(int l_max, double * r, double * sh){
    auto const sph = cart2sph({r[0], r[1], r[2]});
    for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
            sh[(l * l + l + m)] = Ylm(l, m, sph[0], sph[1]);
        }
    }
}
#endif
