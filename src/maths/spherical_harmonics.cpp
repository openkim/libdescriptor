#include "maths/spherical_harmonics.hpp"

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
    double phi = atan2(r_norm[1], r_norm[0]);
    if (phi < 0.0) {
        phi += 2.0 * M_PI;
    }
    return {phi, theta, norm};
}


// *********************************************************************************

// Condon-Shortley phase: (-1)^m
double condon_shortley(int m) {
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

    if (l < 11){
        return Ylmi_analytical(l, m, theta, phi);
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
    // TODO: use pointer based fast sph harm version
    auto const sph = cart2sph({r[0], r[1], r[2]});
    int idx = 0;
    for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
            auto sh_i = (Ylmi(l, m, sph[0], sph[1]));
            sh_real[idx] = sh_i[0];
            sh_imag[idx] = sh_i[1];
            idx++;
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


// *********************************************************************************
// Analytical SH functions *********************************************************
// *********************************************************************************

std::array<double, 2> Y_0_0(double theta, double phi){
    double sph = 1.0/2 * std::sqrt(1/(M_PI) );
    return {sph * 1.0, 0.0};
}


std::array<double, 2> Y_1_neg1(double theta, double phi){
    double sph = 1.0/2 * std::sqrt(3/(2 * M_PI) ) * sin(theta);
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_1_0(double theta, double phi){
    double sph = 1.0/2 * std::sqrt(3/(M_PI) ) * cos(theta);
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_1_1(double theta, double phi){
    double sph = -1.0/2 * std::sqrt(3/(2 * M_PI) ) * sin(theta);
    return {sph * cos(phi), sph * sin(phi)};
}


std::array<double, 2> Y_2_neg2(double theta, double phi){
    double sph = 1.0/4 * std::sqrt(15/(2 * M_PI) ) * std::pow(sin(theta), 2);
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_2_neg1(double theta, double phi){
    double sph = 1.0/2 * std::sqrt(15/(2 * M_PI) ) * sin(theta) * cos(theta);
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_2_0(double theta, double phi){
    double sph = 1.0/4 * std::sqrt(5/(M_PI) ) * (3 * std::pow(cos(theta) ,2) - 1) ;
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_2_1(double theta, double phi){
    double sph = -1.0/2 * std::sqrt(15/(2 * M_PI) ) * sin(theta) * cos(theta);
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_2_2(double theta, double phi){
    double sph = 1.0/4 * std::sqrt(15/(2 * M_PI) ) * std::pow(sin(theta), 2);
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}


std::array<double, 2> Y_3_neg3(double theta, double phi){
    double sph = 1.0/8 * std::sqrt(35/(M_PI) ) * std::pow(sin(theta), 3) ;
    return {sph * cos(-3 * phi), sph * sin(-3 * phi)};
}
std::array<double, 2> Y_3_neg2(double theta, double phi){
    double sph = 1.0/4 * std::sqrt(105/(2 * M_PI) ) * std::pow(sin(theta), 2) * cos(theta) ;
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_3_neg1(double theta, double phi){
    double sph = 1.0/8 * std::sqrt(21/(M_PI) ) * sin(theta) * (5* std::pow(cos(theta) ,2) - 1) ;
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_3_0(double theta, double phi){
    double sph = 1.0/4 * std::sqrt(7/(M_PI) ) * (5* std::pow(cos(theta) ,3) - 3* cos(theta)) ;
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_3_1(double theta, double phi){
    double sph = -1.0/8 * std::sqrt(21/(M_PI) ) * sin(theta) * (5* std::pow(cos(theta) ,2) - 1) ;
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_3_2(double theta, double phi){
    double sph = 1.0/4 * std::sqrt(105/(2 * M_PI) ) * std::pow(sin(theta), 2) * cos(theta) ;
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}
std::array<double, 2> Y_3_3(double theta, double phi){
    double sph = -1.0/8 * std::sqrt(35/(M_PI) ) * std::pow(sin(theta), 3) ;
    return {sph * cos(3 * phi), sph * sin(3 * phi)};
}



std::array<double, 2> Y_4_neg4(double theta, double phi){
    double sph = 3.0/16 * std::sqrt(35/(2 * M_PI) ) * std::pow(sin(theta), 4) ;
    return {sph * cos(-4 * phi), sph * sin(-4 * phi)};
}
std::array<double, 2> Y_4_neg3(double theta, double phi){
    double sph = 3.0/8 * std::sqrt(35/(M_PI) ) * std::pow(sin(theta), 3) * cos(theta) ;
    return {sph * cos(-3 * phi), sph * sin(-3 * phi)};
}
std::array<double, 2> Y_4_neg2(double theta, double phi){
    double sph = 3.0/8 * std::sqrt(5/(2 * M_PI) ) * std::pow(sin(theta), 2) * (7* std::pow(cos(theta) ,2) - 1) ;
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_4_neg1(double theta, double phi){
    double sph = 3.0/8 * std::sqrt(5/(M_PI) ) * sin(theta) * (7* std::pow(cos(theta) ,3) - 3* cos(theta)) ;
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_4_0(double theta, double phi){
    double sph = 3.0/16 * std::sqrt(1/(M_PI) ) * (35* std::pow(cos(theta) ,4) - 30* std::pow(cos(theta) ,2) + 3) ;
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_4_1(double theta, double phi){
    double sph = -3.0/8 * std::sqrt(5/(M_PI) ) * sin(theta) * (7* std::pow(cos(theta) ,3) - 3* cos(theta)) ;
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_4_2(double theta, double phi){
    double sph = 3.0/8 * std::sqrt(5/(2 * M_PI) ) * std::pow(sin(theta), 2) * (7* std::pow(cos(theta) ,2) - 1) ;
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}
std::array<double, 2> Y_4_3(double theta, double phi){
    double sph = -3.0/8 * std::sqrt(35/(M_PI) ) * std::pow(sin(theta), 3) * cos(theta) ;
    return {sph * cos(3 * phi), sph * sin(3 * phi)};
}
std::array<double, 2> Y_4_4(double theta, double phi){
    double sph = 3.0/16 * std::sqrt(35/(2 * M_PI) ) * std::pow(sin(theta), 4) ;
    return {sph * cos(4 * phi), sph * sin(4 * phi)};
}



std::array<double, 2> Y_5_neg5(double theta, double phi){
    double sph = 3.0/32 * std::sqrt(77/(M_PI) ) * std::pow(sin(theta), 5);
    return {sph * cos(-5 * phi), sph * sin(-5 * phi)};
}
std::array<double, 2> Y_5_neg4(double theta, double phi){
    double sph = 3.0/16 * std::sqrt(385/(2 * M_PI) ) * std::pow(sin(theta), 4) * cos(theta);
    return {sph * cos(-4 * phi), sph * sin(-4 * phi)};
}
std::array<double, 2> Y_5_neg3(double theta, double phi){
    double sph = 1.0/32 * std::sqrt(385/(M_PI) ) * std::pow(sin(theta), 3) * (9* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-3 * phi), sph * sin(-3 * phi)};
}
std::array<double, 2> Y_5_neg2(double theta, double phi){
    double sph = 1.0/8 * std::sqrt(1155/(2 * M_PI) ) * std::pow(sin(theta), 2) * (3* std::pow(cos(theta) ,3) - cos(theta));
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_5_neg1(double theta, double phi){
    double sph = 1.0/16 * std::sqrt(165/(2 * M_PI) ) * sin(theta) * (21* std::pow(cos(theta) ,4) - 14* std::pow(cos(theta) ,2) + 1);
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_5_0(double theta, double phi){
    double sph = 1.0/16 * std::sqrt(11/(M_PI) ) * (63* std::pow(cos(theta) ,5) - 70* std::pow(cos(theta) ,3) + 15* cos(theta));
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_5_1(double theta, double phi){
    double sph = -1.0/16 * std::sqrt(165/(2 * M_PI) ) * sin(theta) * (21* std::pow(cos(theta) ,4) - 14* std::pow(cos(theta) ,2) + 1);
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_5_2(double theta, double phi){
    double sph = 1.0/8 * std::sqrt(1155/(2 * M_PI) ) * std::pow(sin(theta), 2) * (3* std::pow(cos(theta) ,3) - cos(theta));
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}
std::array<double, 2> Y_5_3(double theta, double phi){
    double sph = -1.0/32 * std::sqrt(385/(M_PI) ) * std::pow(sin(theta), 3) * (9* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(3 * phi), sph * sin(3 * phi)};
}
std::array<double, 2> Y_5_4(double theta, double phi){
    double sph = 3.0/16 * std::sqrt(385/(2 * M_PI) ) * std::pow(sin(theta), 4) * cos(theta);
    return {sph * cos(4 * phi), sph * sin(4 * phi)};
}
std::array<double, 2> Y_5_5(double theta, double phi){
    double sph = -3.0/32 * std::sqrt(77/(M_PI) ) * std::pow(sin(theta), 5);
    return {sph * cos(5 * phi), sph * sin(5 * phi)};
}



std::array<double, 2> Y_6_neg6(double theta, double phi){
    double sph = 1.0/64 * std::sqrt(3003/(M_PI) ) * std::pow(sin(theta), 6);
    return {sph * cos(-6 * phi), sph * sin(-6 * phi)};
}
std::array<double, 2> Y_6_neg5(double theta, double phi){
    double sph = 3.0/32 * std::sqrt(1001/(M_PI) ) * std::pow(sin(theta), 5) * cos(theta);
    return {sph * cos(-5 * phi), sph * sin(-5 * phi)};
}
std::array<double, 2> Y_6_neg4(double theta, double phi){
    double sph = 3.0/32 * std::sqrt(91/(2 * M_PI) ) * std::pow(sin(theta), 4) * (11* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-4 * phi), sph * sin(-4 * phi)};
}
std::array<double, 2> Y_6_neg3(double theta, double phi){
    double sph = 1.0/32 * std::sqrt(1365/(M_PI) ) * std::pow(sin(theta), 3) * (11* std::pow(cos(theta) ,3) - 3* cos(theta));
    return {sph * cos(-3 * phi), sph * sin(-3 * phi)};
}
std::array<double, 2> Y_6_neg2(double theta, double phi){
    double sph = 1.0/64 * std::sqrt(1365/(M_PI) ) * std::pow(sin(theta), 2) * (33* std::pow(cos(theta) ,4) - 18* std::pow(cos(theta) ,2) + 1);
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_6_neg1(double theta, double phi){
    double sph = 1.0/16 * std::sqrt(273/(2 * M_PI) ) * sin(theta) * (33* std::pow(cos(theta) ,5) - 30* std::pow(cos(theta) ,3) + 5* cos(theta));
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_6_0(double theta, double phi){
    double sph = 1.0/32 * std::sqrt(13/(M_PI) ) * (231* std::pow(cos(theta) ,6) - 315* std::pow(cos(theta) ,4) + 105* std::pow(cos(theta) ,2) - 5);
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_6_1(double theta, double phi){
    double sph = -1.0/16 * std::sqrt(273/(2 * M_PI) ) * sin(theta) * (33* std::pow(cos(theta) ,5) - 30* std::pow(cos(theta) ,3) + 5* cos(theta));
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_6_2(double theta, double phi){
    double sph = 1.0/64 * std::sqrt(1365/(M_PI) ) * std::pow(sin(theta), 2) * (33* std::pow(cos(theta) ,4) - 18* std::pow(cos(theta) ,2) + 1);
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}
std::array<double, 2> Y_6_3(double theta, double phi){
    double sph = -1.0/32 * std::sqrt(1365/(M_PI) ) * std::pow(sin(theta), 3) * (11* std::pow(cos(theta) ,3) - 3* cos(theta));
    return {sph * cos(3 * phi), sph * sin(3 * phi)};
}
std::array<double, 2> Y_6_4(double theta, double phi){
    double sph = 3.0/32 * std::sqrt(91/(2 * M_PI) ) * std::pow(sin(theta), 4) * (11* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(4 * phi), sph * sin(4 * phi)};
}
std::array<double, 2> Y_6_5(double theta, double phi){
    double sph = -3.0/32 * std::sqrt(1001/(M_PI) ) * std::pow(sin(theta), 5) * cos(theta);
    return {sph * cos(5 * phi), sph * sin(5 * phi)};
}
std::array<double, 2> Y_6_6(double theta, double phi){
    double sph = 1.0/64 * std::sqrt(3003/(M_PI) ) * std::pow(sin(theta), 6);
    return {sph * cos(6 * phi), sph * sin(6 * phi)};
}



std::array<double, 2> Y_7_neg7(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(715/(2 * M_PI) ) * std::pow(sin(theta), 7);
    return {sph * cos(-7 * phi), sph * sin(-7 * phi)};
}
std::array<double, 2> Y_7_neg6(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(5005/(M_PI) ) * std::pow(sin(theta), 6) * cos(theta);
    return {sph * cos(-6 * phi), sph * sin(-6 * phi)};
}
std::array<double, 2> Y_7_neg5(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(385/(2 * M_PI) ) * std::pow(sin(theta), 5) * (13* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-5 * phi), sph * sin(-5 * phi)};
}
std::array<double, 2> Y_7_neg4(double theta, double phi){
    double sph = 3.0/32 * std::sqrt(385/(2 * M_PI) ) * std::pow(sin(theta), 4) * (13* std::pow(cos(theta) ,3) - 3* cos(theta));
    return {sph * cos(-4 * phi), sph * sin(-4 * phi)};
}
std::array<double, 2> Y_7_neg3(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(35/(2 * M_PI) ) * std::pow(sin(theta), 3) * (143* std::pow(cos(theta) ,4) - 66* std::pow(cos(theta) ,2) + 3);
    return {sph * cos(-3 * phi), sph * sin(-3 * phi)};
}
std::array<double, 2> Y_7_neg2(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(35/(M_PI) ) * std::pow(sin(theta), 2) * (143* std::pow(cos(theta) ,5) - 110* std::pow(cos(theta) ,3) + 15* cos(theta));
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_7_neg1(double theta, double phi){
    double sph = 1.0/64 * std::sqrt(105/(2 * M_PI) ) * sin(theta) * (429* std::pow(cos(theta) ,6) - 495* std::pow(cos(theta) ,4) + 135* std::pow(cos(theta) ,2) - 5);
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_7_0(double theta, double phi){
    double sph = 1.0/32 * std::sqrt(15/(M_PI) ) * (429* std::pow(cos(theta) ,7) - 693* std::pow(cos(theta) ,5) + 315* std::pow(cos(theta) ,3) - 35* cos(theta));
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_7_1(double theta, double phi){
    double sph = -1.0/64 * std::sqrt(105/(2 * M_PI) ) * sin(theta) * (429* std::pow(cos(theta) ,6) - 495* std::pow(cos(theta) ,4) + 135* std::pow(cos(theta) ,2) - 5);
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_7_2(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(35/(M_PI) ) * std::pow(sin(theta), 2) * (143* std::pow(cos(theta) ,5) - 110* std::pow(cos(theta) ,3) + 15* cos(theta));
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}
std::array<double, 2> Y_7_3(double theta, double phi){
    double sph = -3.0/64 * std::sqrt(35/(2 * M_PI) ) * std::pow(sin(theta), 3) * (143* std::pow(cos(theta) ,4) - 66* std::pow(cos(theta) ,2) + 3);
    return {sph * cos(3 * phi), sph * sin(3 * phi)};
}
std::array<double, 2> Y_7_4(double theta, double phi){
    double sph = 3.0/32 * std::sqrt(385/(2 * M_PI) ) * std::pow(sin(theta), 4) * (13* std::pow(cos(theta) ,3) - 3* cos(theta));
    return {sph * cos(4 * phi), sph * sin(4 * phi)};
}
std::array<double, 2> Y_7_5(double theta, double phi){
    double sph = -3.0/64 * std::sqrt(385/(2 * M_PI) ) * std::pow(sin(theta), 5) * (13* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(5 * phi), sph * sin(5 * phi)};
}
std::array<double, 2> Y_7_6(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(5005/(M_PI) ) * std::pow(sin(theta), 6) * cos(theta);
    return {sph * cos(6 * phi), sph * sin(6 * phi)};
}
std::array<double, 2> Y_7_7(double theta, double phi){
    double sph = -3.0/64 * std::sqrt(715/(2 * M_PI) ) * std::pow(sin(theta), 7);
    return {sph * cos(7 * phi), sph * sin(7 * phi)};
}



std::array<double, 2> Y_8_neg8(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(12155/(2 * M_PI) ) * std::pow(sin(theta), 8);
    return {sph * cos(-8 * phi), sph * sin(-8 * phi)};
}
std::array<double, 2> Y_8_neg7(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(12155/(2 * M_PI) ) * std::pow(sin(theta), 7) * cos(theta);
    return {sph * cos(-7 * phi), sph * sin(-7 * phi)};
}
std::array<double, 2> Y_8_neg6(double theta, double phi){
    double sph = 1.0/128 * std::sqrt(7293/(M_PI) ) * std::pow(sin(theta), 6) * (15* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-6 * phi), sph * sin(-6 * phi)};
}
std::array<double, 2> Y_8_neg5(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(17017/(2 * M_PI) ) * std::pow(sin(theta), 5) * (5* std::pow(cos(theta) ,3) - cos(theta));
    return {sph * cos(-5 * phi), sph * sin(-5 * phi)};
}
std::array<double, 2> Y_8_neg4(double theta, double phi){
    double sph = 3.0/128 * std::sqrt(1309/(2 * M_PI) ) * std::pow(sin(theta), 4) * (65* std::pow(cos(theta) ,4) - 26* std::pow(cos(theta) ,2) + 1);
    return {sph * cos(-4 * phi), sph * sin(-4 * phi)};
}
std::array<double, 2> Y_8_neg3(double theta, double phi){
    double sph = 1.0/64 * std::sqrt(19635/(2 * M_PI) ) * std::pow(sin(theta), 3) * (39* std::pow(cos(theta) ,5) - 26* std::pow(cos(theta) ,3) + 3* cos(theta));
    return {sph * cos(-3 * phi), sph * sin(-3 * phi)};
}
std::array<double, 2> Y_8_neg2(double theta, double phi){
    double sph = 3.0/128 * std::sqrt(595/(M_PI) ) * std::pow(sin(theta), 2) * (143* std::pow(cos(theta) ,6) - 143* std::pow(cos(theta) ,4) + 33* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_8_neg1(double theta, double phi){
    double sph = 3.0/64 * std::sqrt(17/(2 * M_PI) ) * sin(theta) * (715* std::pow(cos(theta) ,7) - 1001* std::pow(cos(theta) ,5) + 385* std::pow(cos(theta) ,3) - 35* cos(theta));
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_8_0(double theta, double phi){
    double sph = 1.0/256 * std::sqrt(17/(M_PI) ) * (6435* std::pow(cos(theta) ,8) - 12012* std::pow(cos(theta) ,6) + 6930* std::pow(cos(theta) ,4) - 1260* std::pow(cos(theta) ,2) + 35);
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_8_1(double theta, double phi){
    double sph = -3.0/64 * std::sqrt(17/(2 * M_PI) ) * sin(theta) * (715* std::pow(cos(theta) ,7) - 1001* std::pow(cos(theta) ,5) + 385* std::pow(cos(theta) ,3) - 35* cos(theta));
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_8_2(double theta, double phi){
    double sph = 3.0/128 * std::sqrt(595/(M_PI) ) * std::pow(sin(theta), 2) * (143* std::pow(cos(theta) ,6) - 143* std::pow(cos(theta) ,4) + 33* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}
std::array<double, 2> Y_8_3(double theta, double phi){
    double sph = -1.0/64 * std::sqrt(19635/(2 * M_PI) ) * std::pow(sin(theta), 3) * (39* std::pow(cos(theta) ,5) - 26* std::pow(cos(theta) ,3) + 3* cos(theta));
    return {sph * cos(3 * phi), sph * sin(3 * phi)};
}
std::array<double, 2> Y_8_4(double theta, double phi){
    double sph = 3.0/128 * std::sqrt(1309/(2 * M_PI) ) * std::pow(sin(theta), 4) * (65* std::pow(cos(theta) ,4) - 26* std::pow(cos(theta) ,2) + 1);
    return {sph * cos(4 * phi), sph * sin(4 * phi)};
}
std::array<double, 2> Y_8_5(double theta, double phi){
    double sph = -3.0/64 * std::sqrt(17017/(2 * M_PI) ) * std::pow(sin(theta), 5) * (5* std::pow(cos(theta) ,3) - cos(theta));
    return {sph * cos(5 * phi), sph * sin(5 * phi)};
}
std::array<double, 2> Y_8_6(double theta, double phi){
    double sph = 1.0/128 * std::sqrt(7293/(M_PI) ) * std::pow(sin(theta), 6) * (15* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(6 * phi), sph * sin(6 * phi)};
}
std::array<double, 2> Y_8_7(double theta, double phi){
    double sph = -3.0/64 * std::sqrt(12155/(2 * M_PI) ) * std::pow(sin(theta), 7) * cos(theta);
    return {sph * cos(7 * phi), sph * sin(7 * phi)};
}
std::array<double, 2> Y_8_8(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(12155/(2 * M_PI) ) * std::pow(sin(theta), 8);
    return {sph * cos(8 * phi), sph * sin(8 * phi)};
}



std::array<double, 2> Y_9_neg9(double theta, double phi){
    double sph = 1.0/512 * std::sqrt(230945/(M_PI) ) * std::pow(sin(theta), 9);
    return {sph * cos(-9 * phi), sph * sin(-9 * phi)};
}
std::array<double, 2> Y_9_neg8(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(230945/(2 * M_PI) ) * std::pow(sin(theta), 8) * cos(theta);
    return {sph * cos(-8 * phi), sph * sin(-8 * phi)};
}
std::array<double, 2> Y_9_neg7(double theta, double phi){
    double sph = 3.0/512 * std::sqrt(13585/(M_PI) ) * std::pow(sin(theta), 7) * (17* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-7 * phi), sph * sin(-7 * phi)};
}
std::array<double, 2> Y_9_neg6(double theta, double phi){
    double sph = 1.0/128 * std::sqrt(40755/(M_PI) ) * std::pow(sin(theta), 6) * (17* std::pow(cos(theta) ,3) - 3* cos(theta));
    return {sph * cos(-6 * phi), sph * sin(-6 * phi)};
}
std::array<double, 2> Y_9_neg5(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(2717/(M_PI) ) * std::pow(sin(theta), 5) * (85* std::pow(cos(theta) ,4) - 30* std::pow(cos(theta) ,2) + 1);
    return {sph * cos(-5 * phi), sph * sin(-5 * phi)};
}
std::array<double, 2> Y_9_neg4(double theta, double phi){
    double sph = 3.0/128 * std::sqrt(95095/(2 * M_PI) ) * std::pow(sin(theta), 4) * (17* std::pow(cos(theta) ,5) - 10* std::pow(cos(theta) ,3) + cos(theta));
    return {sph * cos(-4 * phi), sph * sin(-4 * phi)};
}
std::array<double, 2> Y_9_neg3(double theta, double phi){
    double sph = 1.0/256 * std::sqrt(21945/(M_PI) ) * std::pow(sin(theta), 3) * (221* std::pow(cos(theta) ,6) - 195* std::pow(cos(theta) ,4) + 39* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-3 * phi), sph * sin(-3 * phi)};
}
std::array<double, 2> Y_9_neg2(double theta, double phi){
    double sph = 3.0/128 * std::sqrt(1045/(M_PI) ) * std::pow(sin(theta), 2) * (221* std::pow(cos(theta) ,7) - 273* std::pow(cos(theta) ,5) + 91* std::pow(cos(theta) ,3) - 7* cos(theta));
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_9_neg1(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(95/(2 * M_PI) ) * sin(theta) * (2431* std::pow(cos(theta) ,8) - 4004* std::pow(cos(theta) ,6) + 2002* std::pow(cos(theta) ,4) - 308* std::pow(cos(theta) ,2) + 7);
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_9_0(double theta, double phi){
    double sph = 1.0/256 * std::sqrt(19/(M_PI) ) * (12155* std::pow(cos(theta) ,9) - 25740* std::pow(cos(theta) ,7) + 18018* std::pow(cos(theta) ,5) - 4620* std::pow(cos(theta) ,3) + 315* cos(theta));
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_9_1(double theta, double phi){
    double sph = -3.0/256 * std::sqrt(95/(2 * M_PI) ) * sin(theta) * (2431* std::pow(cos(theta) ,8) - 4004* std::pow(cos(theta) ,6) + 2002* std::pow(cos(theta) ,4) - 308* std::pow(cos(theta) ,2) + 7);
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_9_2(double theta, double phi){
    double sph = 3.0/128 * std::sqrt(1045/(M_PI) ) * std::pow(sin(theta), 2) * (221* std::pow(cos(theta) ,7) - 273* std::pow(cos(theta) ,5) + 91* std::pow(cos(theta) ,3) - 7* cos(theta));
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}
std::array<double, 2> Y_9_3(double theta, double phi){
    double sph = -1.0/256 * std::sqrt(21945/(M_PI) ) * std::pow(sin(theta), 3) * (221* std::pow(cos(theta) ,6) - 195* std::pow(cos(theta) ,4) + 39* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(3 * phi), sph * sin(3 * phi)};
}
std::array<double, 2> Y_9_4(double theta, double phi){
    double sph = 3.0/128 * std::sqrt(95095/(2 * M_PI) ) * std::pow(sin(theta), 4) * (17* std::pow(cos(theta) ,5) - 10* std::pow(cos(theta) ,3) + cos(theta));
    return {sph * cos(4 * phi), sph * sin(4 * phi)};
}
std::array<double, 2> Y_9_5(double theta, double phi){
    double sph = -3.0/256 * std::sqrt(2717/(M_PI) ) * std::pow(sin(theta), 5) * (85* std::pow(cos(theta) ,4) - 30* std::pow(cos(theta) ,2) + 1);
    return {sph * cos(5 * phi), sph * sin(5 * phi)};
}
std::array<double, 2> Y_9_6(double theta, double phi){
    double sph = 1.0/128 * std::sqrt(40755/(M_PI) ) * std::pow(sin(theta), 6) * (17* std::pow(cos(theta) ,3) - 3* cos(theta));
    return {sph * cos(6 * phi), sph * sin(6 * phi)};
}
std::array<double, 2> Y_9_7(double theta, double phi){
    double sph = -3.0/512 * std::sqrt(13585/(M_PI) ) * std::pow(sin(theta), 7) * (17* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(7 * phi), sph * sin(7 * phi)};
}
std::array<double, 2> Y_9_8(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(230945/(2 * M_PI) ) * std::pow(sin(theta), 8) * cos(theta);
    return {sph * cos(8 * phi), sph * sin(8 * phi)};
}
std::array<double, 2> Y_9_9(double theta, double phi){
    double sph = -1.0/512 * std::sqrt(230945/(M_PI) ) * std::pow(sin(theta), 9);
    return {sph * cos(9 * phi), sph * sin(9 * phi)};
}



std::array<double, 2> Y_10_neg10(double theta, double phi){
    double sph = 1.0/1024 * std::sqrt(969969/(M_PI) ) * std::pow(sin(theta), 10);
    return {sph * cos(-10 * phi), sph * sin(-10 * phi)};
}
std::array<double, 2> Y_10_neg9(double theta, double phi){
    double sph = 1.0/512 * std::sqrt(4849845/(M_PI) ) * std::pow(sin(theta), 9) * cos(theta);
    return {sph * cos(-9 * phi), sph * sin(-9 * phi)};
}
std::array<double, 2> Y_10_neg8(double theta, double phi){
    double sph = 1.0/512 * std::sqrt(255255/(2 * M_PI) ) * std::pow(sin(theta), 8) * (19* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-8 * phi), sph * sin(-8 * phi)};
}
std::array<double, 2> Y_10_neg7(double theta, double phi){
    double sph = 3.0/512 * std::sqrt(85085/(M_PI) ) * std::pow(sin(theta), 7) * (19* std::pow(cos(theta) ,3) - 3* cos(theta));
    return {sph * cos(-7 * phi), sph * sin(-7 * phi)};
}
std::array<double, 2> Y_10_neg6(double theta, double phi){
    double sph = 3.0/1024 * std::sqrt(5005/(M_PI) ) * std::pow(sin(theta), 6) * (323* std::pow(cos(theta) ,4) - 102* std::pow(cos(theta) ,2) + 3);
    return {sph * cos(-6 * phi), sph * sin(-6 * phi)};
}
std::array<double, 2> Y_10_neg5(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(1001/(M_PI) ) * std::pow(sin(theta), 5) * (323* std::pow(cos(theta) ,5) - 170* std::pow(cos(theta) ,3) + 15* cos(theta));
    return {sph * cos(-5 * phi), sph * sin(-5 * phi)};
}
std::array<double, 2> Y_10_neg4(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(5005/(2 * M_PI) ) * std::pow(sin(theta), 4) * (323* std::pow(cos(theta) ,6) - 255* std::pow(cos(theta) ,4) + 45* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(-4 * phi), sph * sin(-4 * phi)};
}
std::array<double, 2> Y_10_neg3(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(5005/(M_PI) ) * std::pow(sin(theta), 3) * (323* std::pow(cos(theta) ,7) - 357* std::pow(cos(theta) ,5) + 105* std::pow(cos(theta) ,3) - 7* cos(theta));
    return {sph * cos(-3 * phi), sph * sin(-3 * phi)};
}
std::array<double, 2> Y_10_neg2(double theta, double phi){
    double sph = 3.0/512 * std::sqrt(385/(2 * M_PI) ) * std::pow(sin(theta), 2) * (4199* std::pow(cos(theta) ,8) - 6188* std::pow(cos(theta) ,6) + 2730* std::pow(cos(theta) ,4) - 364* std::pow(cos(theta) ,2) + 7);
    return {sph * cos(-2 * phi), sph * sin(-2 * phi)};
}
std::array<double, 2> Y_10_neg1(double theta, double phi){
    double sph = 1.0/256 * std::sqrt(1155/(2 * M_PI) ) * sin(theta) * (4199* std::pow(cos(theta) ,9) - 7956* std::pow(cos(theta) ,7) + 4914* std::pow(cos(theta) ,5) - 1092* std::pow(cos(theta) ,3) + 63* cos(theta));
    return {sph * cos(-phi), sph * sin(-phi)};
}
std::array<double, 2> Y_10_0(double theta, double phi){
    double sph = 1.0/512 * std::sqrt(21/(M_PI) ) * (46189 * std::pow(cos(theta), 10) - 109395* std::pow(cos(theta) ,8) + 90090* std::pow(cos(theta) ,6) - 30030* std::pow(cos(theta) ,4) + 3465* std::pow(cos(theta) ,2) - 63);
    return {sph * 1.0, 0.0};
}
std::array<double, 2> Y_10_1(double theta, double phi){
    double sph = -1.0/256 * std::sqrt(1155/(2 * M_PI) ) * sin(theta) * (4199* std::pow(cos(theta) ,9) - 7956* std::pow(cos(theta) ,7) + 4914* std::pow(cos(theta) ,5) - 1092* std::pow(cos(theta) ,3) + 63* cos(theta));
    return {sph * cos(phi), sph * sin(phi)};
}
std::array<double, 2> Y_10_2(double theta, double phi){
    double sph = 3.0/512 * std::sqrt(385/(2 * M_PI) ) * std::pow(sin(theta), 2) * (4199* std::pow(cos(theta) ,8) - 6188* std::pow(cos(theta) ,6) + 2730* std::pow(cos(theta) ,4) - 364* std::pow(cos(theta) ,2) + 7);
    return {sph * cos(2 * phi), sph * sin(2 * phi)};
}
std::array<double, 2> Y_10_3(double theta, double phi){
    double sph = -3.0/256 * std::sqrt(5005/(M_PI) ) * std::pow(sin(theta), 3) * (323* std::pow(cos(theta) ,7) - 357* std::pow(cos(theta) ,5) + 105* std::pow(cos(theta) ,3) - 7* cos(theta));
    return {sph * cos(3 * phi), sph * sin(3 * phi)};
}
std::array<double, 2> Y_10_4(double theta, double phi){
    double sph = 3.0/256 * std::sqrt(5005/(2 * M_PI) ) * std::pow(sin(theta), 4) * (323* std::pow(cos(theta) ,6) - 255* std::pow(cos(theta) ,4) + 45* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(4 * phi), sph * sin(4 * phi)};
}
std::array<double, 2> Y_10_5(double theta, double phi){
    double sph = -3.0/256 * std::sqrt(1001/(M_PI) ) * std::pow(sin(theta), 5) * (323* std::pow(cos(theta) ,5) - 170* std::pow(cos(theta) ,3) + 15* cos(theta));
    return {sph * cos(5 * phi), sph * sin(5 * phi)};
}
std::array<double, 2> Y_10_6(double theta, double phi){
    double sph = 3.0/1024 * std::sqrt(5005/(M_PI) ) * std::pow(sin(theta), 6) * (323* std::pow(cos(theta) ,4) - 102* std::pow(cos(theta) ,2) + 3);
    return {sph * cos(6 * phi), sph * sin(6 * phi)};
}
std::array<double, 2> Y_10_7(double theta, double phi){
    double sph = -3.0/512 * std::sqrt(85085/(M_PI) ) * std::pow(sin(theta), 7) * (19* std::pow(cos(theta) ,3) - 3* cos(theta));
    return {sph * cos(7 * phi), sph * sin(7 * phi)};
}
std::array<double, 2> Y_10_8(double theta, double phi){
    double sph = 1.0/512 * std::sqrt(255255/(2 * M_PI) ) * std::pow(sin(theta), 8) * (19* std::pow(cos(theta) ,2) - 1);
    return {sph * cos(8 * phi), sph * sin(8 * phi)};
}
std::array<double, 2> Y_10_9(double theta, double phi){
    double sph = -1.0/512 * std::sqrt(4849845/(M_PI) ) * std::pow(sin(theta), 9) * cos(theta);
    return {sph * cos(9 * phi), sph * sin(9 * phi)};
}
std::array<double, 2> Y_10_10(double theta, double phi){
    double sph = 1.0/1024 * std::sqrt(969969/(M_PI) ) * std::pow(sin(theta), 10);
    return {sph * cos(10 * phi), sph * sin(10 * phi)};
}


std::array<double, 2> Ylmi_analytical(int l,int m, double theta, double phi){
    switch (l) {
        case 0: {
            return Y_0_0(theta, phi);
        }
        case 1: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -1: return Y_1_neg1(theta, phi);
                case 0: return Y_1_0(theta, phi);
                case 1: return Y_1_1(theta, phi);
            }
        }
        case 2: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -2: return Y_2_neg2(theta, phi);
                case -1: return Y_2_neg1(theta, phi);
                case 0: return Y_2_0(theta, phi);
                case 1: return Y_2_1(theta, phi);
                case 2: return Y_2_2(theta, phi);
            }
        }
        case 3: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -3:return Y_3_neg3(theta, phi);
                case -2:return Y_3_neg2(theta, phi);
                case -1:return Y_3_neg1(theta, phi);
                case 0:return Y_3_0(theta, phi);
                case 1:return Y_3_1(theta, phi);
                case 2:return Y_3_2(theta, phi);
                case 3:return Y_3_3(theta, phi);
            }
        }
        case 4: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -4: return Y_4_neg4(theta, phi);
                case -3: return Y_4_neg3(theta, phi);
                case -2: return Y_4_neg2(theta, phi);
                case -1: return Y_4_neg1(theta, phi);
                case 0: return Y_4_0(theta, phi);
                case 1: return Y_4_1(theta, phi);
                case 2: return Y_4_2(theta, phi);
                case 3: return Y_4_3(theta, phi);
                case 4: return Y_4_4(theta, phi);
            }
        }
        case 5: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -5: return Y_5_neg5(theta, phi);
                case -4: return Y_5_neg4(theta, phi);
                case -3: return Y_5_neg3(theta, phi);
                case -2: return Y_5_neg2(theta, phi);
                case -1: return Y_5_neg1(theta, phi);
                case 0: return Y_5_0(theta, phi);
                case 1: return Y_5_1(theta, phi);
                case 2: return Y_5_2(theta, phi);
                case 3: return Y_5_3(theta, phi);
                case 4: return Y_5_4(theta, phi);
                case 5: return Y_5_5(theta, phi);
            }
        }
        case 6: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -6: return Y_6_neg6(theta, phi);
                case -5: return Y_6_neg5(theta, phi);
                case -4: return Y_6_neg4(theta, phi);
                case -3: return Y_6_neg3(theta, phi);
                case -2: return Y_6_neg2(theta, phi);
                case -1: return Y_6_neg1(theta, phi);
                case 0: return Y_6_0(theta, phi);
                case 1: return Y_6_1(theta, phi);
                case 2: return Y_6_2(theta, phi);
                case 3: return Y_6_3(theta, phi);
                case 4: return Y_6_4(theta, phi);
                case 5: return Y_6_5(theta, phi);
                case 6: return Y_6_6(theta, phi);
            }
        }
        case 7: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -7: return Y_7_neg7(theta, phi);
                case -6: return Y_7_neg6(theta, phi);
                case -5: return Y_7_neg5(theta, phi);
                case -4: return Y_7_neg4(theta, phi);
                case -3: return Y_7_neg3(theta, phi);
                case -2: return Y_7_neg2(theta, phi);
                case -1: return Y_7_neg1(theta, phi);
                case 0: return Y_7_0(theta, phi);
                case 1: return Y_7_1(theta, phi);
                case 2: return Y_7_2(theta, phi);
                case 3: return Y_7_3(theta, phi);
                case 4: return Y_7_4(theta, phi);
                case 5: return Y_7_5(theta, phi);
                case 6: return Y_7_6(theta, phi);
                case 7: return Y_7_7(theta, phi);
            }
        }
        case 8: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -8: return Y_8_neg8(theta, phi);
                case -7: return Y_8_neg7(theta, phi);
                case -6: return Y_8_neg6(theta, phi);
                case -5: return Y_8_neg5(theta, phi);
                case -4: return Y_8_neg4(theta, phi);
                case -3: return Y_8_neg3(theta, phi);
                case -2: return Y_8_neg2(theta, phi);
                case -1: return Y_8_neg1(theta, phi);
                case 0: return Y_8_0(theta, phi);
                case 1: return Y_8_1(theta, phi);
                case 2: return Y_8_2(theta, phi);
                case 3: return Y_8_3(theta, phi);
                case 4: return Y_8_4(theta, phi);
                case 5: return Y_8_5(theta, phi);
                case 6: return Y_8_6(theta, phi);
                case 7: return Y_8_7(theta, phi);
                case 8: return Y_8_8(theta, phi);
            }
        }
        case 9: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -9: return Y_9_neg9(theta, phi);
                case -8: return Y_9_neg8(theta, phi);
                case -7: return Y_9_neg7(theta, phi);
                case -6: return Y_9_neg6(theta, phi);
                case -5: return Y_9_neg5(theta, phi);
                case -4: return Y_9_neg4(theta, phi);
                case -3: return Y_9_neg3(theta, phi);
                case -2: return Y_9_neg2(theta, phi);
                case -1: return Y_9_neg1(theta, phi);
                case 0: return Y_9_0(theta, phi);
                case 1: return Y_9_1(theta, phi);
                case 2: return Y_9_2(theta, phi);
                case 3: return Y_9_3(theta, phi);
                case 4: return Y_9_4(theta, phi);
                case 5: return Y_9_5(theta, phi);
                case 6: return Y_9_6(theta, phi);
                case 7: return Y_9_7(theta, phi);
                case 8: return Y_9_8(theta, phi);
                case 9: return Y_9_9(theta, phi);
            }
        }
        case 10: {
            switch(m){
                default: throw std::runtime_error("Invalid value of m.");
                case -10: return Y_10_neg10(theta, phi);
                case -9: return Y_10_neg9(theta, phi);
                case -8: return Y_10_neg8(theta, phi);
                case -7: return Y_10_neg7(theta, phi);
                case -6: return Y_10_neg6(theta, phi);
                case -5: return Y_10_neg5(theta, phi);
                case -4: return Y_10_neg4(theta, phi);
                case -3: return Y_10_neg3(theta, phi);
                case -2: return Y_10_neg2(theta, phi);
                case -1: return Y_10_neg1(theta, phi);
                case 0: return Y_10_0(theta, phi);
                case 1: return Y_10_1(theta, phi);
                case 2: return Y_10_2(theta, phi);
                case 3: return Y_10_3(theta, phi);
                case 4: return Y_10_4(theta, phi);
                case 5: return Y_10_5(theta, phi);
                case 6: return Y_10_6(theta, phi);
                case 7: return Y_10_7(theta, phi);
                case 8: return Y_10_8(theta, phi);
                case 9: return Y_10_9(theta, phi);
                case 10: return Y_10_10(theta, phi);
            }
        }
        default:{
            throw std::runtime_error("Only l <= 10 have analytical expressions.");
        }
    }
}