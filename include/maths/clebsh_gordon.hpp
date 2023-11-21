
#ifndef CLEBSH_GORDON_HPP
#define CLEBSH_GORDON_HPP

#define _USE_MATH_DEFINES
// Clebsh-Gordan coefficients *******************************************************
// Weigner 3j symbols
// Based on https://earlelab.rit.albany.edu/documents/Lai369j.pdf and ETAM ME Rose
// using racah formula, for some reason weigner formula fails for some cases (usually m > 5)
// lgamman function suggested by: Kobi (2023). Wigner3j symbol (https://www.mathworks.com/matlabcentral/fileexchange/20619-wigner3j-symbol), MATLAB Central File Exchange. Retrieved July 19, 2023.
double weigner_3j(double j1, double j2, double j3, double m1, double m2, double m3);

// Clebsh-Gordan coefficients, usage <j1m1,j2m2|j3m3>
double clebsh_gordon(double j1, double j2, double j3, double m1, double m2, double m3);

#endif //CLEBSH_GORDON_HPP