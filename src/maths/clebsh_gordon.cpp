#include <cmath>
#include "maths/gamma.hpp"
#include "maths/clebsh_gordon.hpp"

#define SIGN_POWER(m) ((m) % 2 == 0 ? 1 : -1)

// Clebsh-Gordan coefficients *******************************************************
// Weigner 3j symbols
// Based on https://earlelab.rit.albany.edu/documents/Lai369j.pdf and ETAM ME Rose
// using racah formula, for some reason weigner formula fails for some cases (usually m > 5)
// lgamman function suggested by: Kobi (2023). Wigner3j symbol (https://www.mathworks.com/matlabcentral/fileexchange/20619-wigner3j-symbol), MATLAB Central File Exchange. Retrieved July 19, 2023.
double weigner_3j(double j1, double j2, double j3, double m1, double m2, double m3){
    // check if the 3j symbol is trivially zero
    double w3j = 0.0;
    if (m1 + m2 + m3 != 0 ||
        j3 < std::abs(j1 - j2) || j3 > j1 + j2 ||
        std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m3) > j3 ||
        j1 < 0 || j2 < 0 || j3 < 0 ||
        (m1 == 0 && m2 == 0 && m3 == 0 && ((static_cast<int>(j1 + j2 + j3) % 2) != 0))) {
        w3j = 0.0;
    } else {
        double t1 = j2 - m1 - j3;
        double t2 = j1 + m2 - j3;
        double t3 = j1 + j2 - j3;
        double t4 = j1 - m1;
        double t5 = j2 + m2;
        int tmin = static_cast<int>(std::max(0.0, std::max(t1, t2)));
        int tmax = static_cast<int>(std::min(t3, std::min(t4, t5)));
        long double s = 0; // keeping long double is needed for high precision
        // explore using doule to ensure smooth AD (and maybe faster?)
        for (int t = tmin; t <= tmax; t++) {
            s += SIGN_POWER(t) /
                 (std::exp(lgamma(1 + t) + lgamma(1 + t - t1) + lgamma(1 + t - t2) + lgamma(1 + t3 - t) +
                      lgamma(1 + t4 - t) +
                      lgamma(1 + t5 - t)));
        }
        s *= SIGN_POWER(static_cast<int>(j1 - j2 - m3)) *
             sqrt(std::exp(lgamma(1 + j1 + j2 - j3) + lgamma(1 + j1 - j2 + j3) +
                      lgamma(1 + -j1 + j2 + j3) - lgamma(1 + j1 + j2 + j3 + 1) + lgamma(1 + j1 + m1) +
                      lgamma(1 + j1 - m1) + lgamma(1 + j2 + m2) +
                      lgamma(1 + j2 - m2) + lgamma(1 + j3 + m3) + lgamma(1 + j3 - m3)));
        w3j = static_cast<double>(s);
    }
    return w3j;
}

// Clebsh-Gordan coefficients, usage <j1m1,j2m2|j3m3>
double clebsh_gordon(double j1, double j2, double j3, double m1, double m2, double m3){
    return SIGN_POWER(static_cast<int>(j1-j2-m3)) * std::sqrt(2*j3+1) * weigner_3j(j1, j2, j3, m1, m2, -m3);
}
