#include <limits>
#include <cmath>

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
            -1.142022680371168e0, 6.5165112670737e-3,3.087090173086e-4, -3.4706269649e-6,
            6.9437664e-9, 3.67795e-11, -1.356e-13};
    static double c2[] = {
            1.843740587300905e0, -7.68528408447867e-2,1.2719271366546e-3, -4.9717367042e-6,
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
    return (rimu * ril1) / ril;
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

