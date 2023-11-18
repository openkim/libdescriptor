#ifndef BESSEL_FUNCTIONS_HPP
#define BESSEL_FUNCTIONS_HPP
#define _USE_MATH_DEFINES

#include <limits>
#include <cmath>
#include <stdexcept>

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

inline double chebev(const double *c, const int m, const double x) {
    double d = 0.0, dd = 0.0, sv;
    int j;
    for (j = m - 1; j > 0; j--) {
        sv = d;
        d = 2. * x * d - dd + c[j];
        dd = sv;
    }
    return x * d - dd + 0.5 * c[0];
}

// From numerical recipies in C++ 3rd edition
// omiting the besselk and derivative parts
double bessel_I(const double nu, const double x) {
    const int MAXIT = 10000;
    const double EPS = std::numeric_limits<double>::epsilon();
    const double FPMIN = std::numeric_limits<double>::min() / EPS;
    const double XMIN = 2.0;
    double a, a1, b, c, d, del, del1, delh, dels, e, f, fact, fact2, ff, gam1, gam2, gammi, gampl, h, p, pimu, q, q1,
            q2, qnew, ril, ril1, rimu, rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2, xi, xi2, xmu, xmu2, xx;
    int i, l, nl;

    static double c1[] = {
            -1.142022680371168e0, 6.5165112670737e-3, 3.087090173086e-4, -3.4706269649e-6,
            6.9437664e-9, 3.67795e-11, -1.356e-13};
    static double c2[] = {
            1.843740587300905e0, -7.68528408447867e-2, 1.2719271366546e-3, -4.9717367042e-6,
            -3.31261198e-8, 2.423096e-10, -1.702e-13, -1.49e-15};

    if (x <= 0.0 || nu < 0.0) {
        throw std::runtime_error("bad arguments in besselik");
    }
    // nl = static_cast<int>(nu + 0.5);
    nl = static_cast<int>(std::lround(nu + 0.5));
    xmu = nu - nl;
    xmu2 = xmu * xmu;
    xi = 1.0 / x;
    xi2 = 2.0 * xi;
    h = nu * xi;
    if (h < FPMIN) h = FPMIN;
    b = xi2 * nu;
    d = 0.0;
    c = h;
    for (i = 0; i < MAXIT; i++) {
        b += xi2;
        d = 1.0 / (b + d);
        c = b + 1.0 / c;
        del = c * d;
        h = del * h;
        if (std::abs(del - 1.0) <= EPS) break;
    }
    if (i >= MAXIT) {
        throw std::runtime_error("x too large in besselik; try asymptotic expansion");
    }
    ril = FPMIN;
    ripl = h * ril;
    ril1 = ril;
    rip1 = ripl;
    fact = nu * xi;
    for (l = nl - 1; l >= 0; l--) {
        ritemp = fact * ril + ripl;
        fact -= xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
    }
    f = ripl / ril;
    if (x < XMIN) {
        x2 = 0.5 * x;
        pimu = M_PI * xmu;
        fact = (std::abs(pimu) < EPS ? 1.0 : pimu / std::sin(pimu));
        d = -std::log(x2);
        e = xmu * d;
        fact2 = (std::abs(e) < EPS ? 1.0 : std::sinh(e) / e);
        xx = 8.0 * (xmu * xmu) - 1.0;
        gam1 = chebev(c1, 7, xx); // double precision run summation for 7 and 8
        gam2 = chebev(c2, 8, xx);
        gampl = gam2 - xmu * gam1;
        gammi = gam2 + xmu * gam1;
        ff = fact * (gam1 * std::cosh(e) + gam2 * fact2 * d);
        sum = ff;
        e = std::exp(e);
        p = 0.5 * e / gampl;
        q = 0.5 / (e * gammi);
        c = 1.0;
        d = x2 * x2;
        sum1 = p;
        for (i = 1; i <= MAXIT; i++) {
            ff = (i * ff + p + q) / (i * i - xmu2);
            c *= (d / i);
            p /= (i - xmu);
            q /= (i + xmu);
            del = c * ff;
            sum += del;
            del1 = c * (p - i * ff);
            sum1 += del1;
            if (std::abs(del) < std::abs(sum) * EPS) {
                break;
            }
        }
        if (i > MAXIT) throw std::runtime_error("bessk series failed to converge");
        rkmu = sum;
        rk1 = sum1 * xi2;
    } else {
        b = 2.0 * (1.0 + x);
        d = 1.0 / b;
        h = delh = d;
        q1 = 0.0;
        q2 = 1.0;
        a1 = 0.25 - xmu2;
        q = c = a1;
        a = -a1;
        s = 1.0 + q * delh;
        for (i = 1; i < MAXIT; i++) {
            a -= 2 * i;
            c = -a * c / (i + 1.0);
            qnew = (q1 - b * q2) / a;
            q1 = q2;
            q2 = qnew;
            q += c * qnew;
            b += 2.0;
            d = 1.0 / (b + a * d);
            delh = (b * d - 1.0) * delh;
            h += delh;
            dels = q * delh;
            s += dels;
            if (std::abs(dels / s) <= EPS) {
                break;
            }
        }
        if (i >= MAXIT) {
            throw std::runtime_error("besselik: failure to converge in cf2");
        }
        h = a1 * h;
        rkmu = std::sqrt(M_PI / (2.0 * x)) * std::exp(-x) / s;
        rk1 = rkmu * (xmu + x + 0.5 - h) * xi;
    }
    rkmup = xmu * xi * rkmu - rk1;
    rimu = xi / (f * rkmu - rkmup);
    return (rimu * ril1) / ril; //po
    // double ipo=(rimu*rip1)/ril;
    // for (i=1;i <= nl;i++) {
    // rktemp=(xmu+i)*xi2*rk1+rkmu;
    // rkmu=rk1;
    // rk1=rktemp;
    // }
    // double ko=rkmu;
    // double kpo=nu*xi*rkmu-rk1;
    // double xik = x;
    // double nuik = nu;
}


double bessel_J(const double nu, const double x) {
    const int MAXIT = 10000;
    const double EPS = std::numeric_limits<double>::epsilon();
    const double FPMIN = std::numeric_limits<double>::min() / EPS;
    const double XMIN = 2.0;
    double a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2, fact3, ff, gam, gam1, gam2,
    gammi, gampl, h, p, pimu, pimu2, q, r, rjl, rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
    temp, w, x2, xi, xi2, xmu, xmu2, xx;
    int i, isign, l, nl;
    const double c1[7] = {-1.142022680371168e0, 6.5165112670737e-3,
                          3.087090173086e-4, -3.4706269649e-6, 6.9437664e-9, 3.67795e-11,
                          -1.356e-13};
    const double c2[8] = {1.843740587300905e0, -7.68528408447867e-2,
                          1.2719271366546e-3, -4.9717367042e-6, -3.31261198e-8, 2.423096e-10,
                          -1.702e-13, -1.49e-15};

    if (x <= 0.0 || nu < 0.0) throw std::runtime_error("bad arguments in besseljy");
    nl = (x < XMIN ? static_cast<int>(std::lround(nu + 0.5)) : std::max(0, int(nu - x + 1.5)));
    xmu = nu - nl;
    xmu2 = xmu * xmu;
    xi = 1.0 / x;
    xi2 = 2.0 * xi;
    w = xi2 / M_PI;
    isign = 1;
    h = nu * xi;
    if (h < FPMIN) h = FPMIN;
    b = xi2 * nu;
    d = 0.0;
    c = h;
    for (i = 0; i < MAXIT; i++) {
        b += xi2;
        d = b - d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = b - 1.0 / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = c * d;
        h = del * h;
        if (d < 0.0) isign = -isign;
        if (std::abs(del - 1.0) <= EPS) break;
    }
    if (i >= MAXIT)
        throw std::runtime_error("x too large in besseljy; try asymptotic expansion");
    rjl = isign * FPMIN;
    rjpl = h * rjl;
    rjl1 = rjl;
    rjp1 = rjpl;
    fact = nu * xi;
    for (l = nl - 1; l >= 0; l--) {
        rjtemp = fact * rjl + rjpl;
        fact -= xi;
        rjpl = fact * rjtemp - rjl;
        rjl = rjtemp;
    }
    if (rjl == 0.0) rjl = EPS;
    f = rjpl / rjl;
    if (x < XMIN) {
        x2 = 0.5 * x;
        pimu = M_PI * xmu;
        fact = (std::abs(pimu) < EPS ? 1.0 : pimu / std::sin(pimu));
        d = -std::log(x2);
        e = xmu * d;
        fact2 = (std::abs(e) < EPS ? 1.0 : sinh(e) / e);
        xx = 8.0 * (xmu * xmu) - 1.0;
        gam1 = chebev(c1, 7, xx);
        gam2 = chebev(c2, 8, xx);
        gampl = gam2 - xmu * gam1;
        gammi = gam2 + xmu * gam1;
        ff = 2.0 / M_PI * fact * (gam1 * std::cosh(e) + gam2 * fact2 * d);
        e = std::exp(e);
        p = e / (gampl * M_PI);
        q = 1.0 / (e * M_PI * gammi);
        pimu2 = 0.5 * pimu;
        fact3 = (std::abs(pimu2) < EPS ? 1.0 : std::sin(pimu2) / pimu2);
        r = M_PI * pimu2 * fact3 * fact3;
        c = 1.0;
        d = -x2 * x2;
        sum = ff + r * q;
        sum1 = p;
        for (i = 1; i <= MAXIT; i++) {
            ff = (i * ff + p + q) / (i * i - xmu2);
            c *= (d / i);
            p /= (i - xmu);
            q /= (i + xmu);
            del = c * (ff + r * q);
            sum += del;
            del1 = c * p - i * del;
            sum1 += del1;
            if (std::abs(del) < (1.0 + std::abs(sum)) * EPS) break;
        }
        if (i > MAXIT) throw std::runtime_error("bessy series failed to converge");
        rymu = -sum;
        ry1 = -sum1 * xi2;
        rymup = xmu * xi * rymu - ry1;
        rjmu = w / (rymup - f * rymu);
    } else {
        a = 0.25 - xmu2;
        p = -0.5 * xi;
        q = 1.0;
        br = 2.0 * x;
        bi = 2.0;
        fact = a * xi / (p * p + q * q);
        cr = br + q * fact;
        ci = bi + p * fact;
        den = br * br + bi * bi;
        dr = br / den;
        di = -bi / den;
        dlr = cr * dr - ci * di;
        dli = cr * di + ci * dr;
        temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        for (i = 1; i < MAXIT; i++) {
            a += 2 * i;
            bi += 2.0;
            dr = a * dr + br;
            di = a * di + bi;
            if (std::abs(dr) + std::abs(di) < FPMIN) dr = FPMIN;
            fact = a / (cr * cr + ci * ci);
            cr = br + cr * fact;
            ci = bi - ci * fact;
            if (std::abs(cr) + std::abs(ci) < FPMIN) cr = FPMIN;
            den = dr * dr + di * di;
            dr /= den;
            di /= -den;
            dlr = cr * dr - ci * di;
            dli = cr * di + ci * dr;
            temp = p * dlr - q * dli;
            q = p * dli + q * dlr;
            p = temp;
            if (std::abs(dlr - 1.0) + std::abs(dli) <= EPS) break;
        }
        if (i >= MAXIT) throw std::runtime_error("cf2 failed in besseljy");
        gam = (p - f) / q;
        rjmu = std::sqrt(w / ((p - f) * gam + q));
        rjmu = SIGN(rjmu, rjl);
        rymu = rjmu * gam;
        rymup = rymu * (p + q / gam);
        // ry1 = xmu * xi * rymu - rymup;
    }
    fact = rjmu / rjl;
    return rjl1 * fact;
    // jo=rjl1*fact;
    // jpo=rjp1*fact;
    // for (i=1;i<=nl;i++) {
    // rytemp=(xmu+i)*xi2*ry1-rymu;
    // rymu=ry1;
    // ry1=rytemp;
}
// yo=rymu;
// ypo=nu*xi*rymu-ry1;
// xjy = x;
// nujy = nu;

double spherical_in(double n, double x) {
    //Modified spherical Bessel function of the first kind
    return bessel_I(n + 0.5, x) * sqrt(M_PI / (2 * x));
}
double spherical_in(int n, double x){
    return spherical_in(static_cast<double>(n), x);
}

double spherical_jn(double n, double x) {
    //Spherical Bessel function of the first kind
    return bessel_J(n + 0.5, x) * sqrt(M_PI / (2 * x));
}
double spherical_jn(int n, double x){
    return spherical_jn(static_cast<double>(n), x);
}

double halleys_root(double l, double lwr_bnd, double upr_bnd){
    double TOL1=2.2204e-13, TOL2=2.2204e-14;
    int MAXIT=100000;

    double x=(lwr_bnd + upr_bnd)/2.;
    double l1 = l+1;

    int i = 0;
    while (i < MAXIT){
        double a = spherical_jn(l, x);
        double b = spherical_jn(l1, x);

        if (std::abs(a) < TOL2) {
            break;
        }

        double x2 = x*x;
        double dx = -2*x*a*(l*a-x*b) / ((l*l1+x2)*a*a - 2*x*(2*l+1)*a*b + 2*x2*b*b);
        if (std::abs(dx) < TOL1) {
            break;
        }

        if (dx > 0) {
            lwr_bnd = x;
        } else {
            upr_bnd = x;
        }

        if ((upr_bnd-x) < dx || dx < (lwr_bnd-x)) {
            dx = (upr_bnd - lwr_bnd) / 2. - x;
        }

        x=x+dx;
        i++;
    }

    if (i>MAXIT-1) {
        throw std::runtime_error("halleys_root() failed to converge.");
    }
    return x;
}

double halleys_root(int l, double lwr_bnd, double upr_bnd){
    return halleys_root(static_cast<double>(l), lwr_bnd, upr_bnd);
}


void spherical_jn_zeros(int n_max, double * u_all ){
    // Expected u_all to be a (n_max+2, n_max+1)
    for (int l = 0; l < n_max + 2; ++l) {
        u_all[l * (n_max + 1)] = M_PI * (l+1);
    } 

    // call Halley's Method
    for (int l = 1; l < n_max + 1; l++) {
        for (int n = 0; n < n_max - l + 2; n++) {
            u_all[n * (n_max + 1) + l] = halleys_root(static_cast<double>(l), 
            u_all[n * (n_max + 1) + (l  - 1)], 
            u_all[(n + 1) * (n_max + 1) + (l - 1)]);
        }
    }
}

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