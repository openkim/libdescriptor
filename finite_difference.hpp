// Finite Diff from Boost
// Removed complex number bits, and other boost namespaces
// Now it is simple header only 
//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_FINITE_DIFFERENCE_HPP
#define BOOST_FINITE_DIFFERENCE_HPP

#include <cmath>
#include <limits>

//typedef void (*f_ptr) (double *, double *);
/*
 * Performs numerical differentiation by finite-differences.
 *
 * All numerical differentiation using finite-differences are ill-conditioned, and these routines are no exception.
 * A simple argument demonstrates that the error is unbounded as h->0.
 * Take the one sides finite difference formula f'(x) = (f(x+h)-f(x))/h.
 * The evaluation of f induces an error as well as the error from the finite-difference approximation, giving
 * |f'(x) - (f(x+h) -f(x))/h| < h|f''(x)|/2 + (|f(x)|+|f(x+h)|)eps/h =: g(h), where eps is the unit roundoff for the type.
 * It is reasonable to choose h in a way that minimizes the maximum error bound g(h).
 * The value of h that minimizes g is h = sqrt(2eps(|f(x)| + |f(x+h)|)/|f''(x)|), and for this value of h the error bound is
 * sqrt(2eps(|f(x+h) +f(x)||f''(x)|)).
 * In fact it is not necessary to compute the ratio (|f(x+h)| + |f(x)|)/|f''(x)|; the error bound of ~\sqrt{\epsilon} still holds if we set it to one.
 *
 *
 * For more details on this method of analysis, see
 *
 * http://www.uio.no/studier/emner/matnat/math/MAT-INF1100/h08/kompendiet/diffint.pdf
 * http://web.archive.org/web/20150420195907/http://www.uio.no/studier/emner/matnat/math/MAT-INF1100/h08/kompendiet/diffint.pdf
 *
 *
 * It can be shown on general grounds that when choosing the optimal h, the maximum error in f'(x) is ~(|f(x)|eps)^k/k+1|f^(k-1)(x)|^1/k+1.
 * From this we can see that full precision can be recovered in the limit k->infinity.
 *
 * References:
 *
 * 1) Fornberg, Bengt. "Generation of finite difference formulas on arbitrarily spaced grids." Mathematics of computation 51.184 (1988): 699-706.
 */

namespace numdiff {

    double make_xph_representable(double x, double h) {
        using std::numeric_limits;
        // Redefine h so that x + h is representable. Not using this trick leads to large error.
        // The compiler flag -ffast-math evaporates these operations . . .
        double temp = x + h;
        h = temp - x;
        // Handle the case x + h == x:
        if (h == 0) {
            h = std::nextafter(x, (numeric_limits<double>::max)()) - x;
        }
        return h;
    }

    template<class F>
    void vec_finite_difference_derivative(F& f, double *x, int pos, int input_size, int output_size, double *y) {
        using std::sqrt;
        using std::pow;
        using std::abs;
        using std::numeric_limits;

        const double eps = (numeric_limits<double>::epsilon)();
        double h = pow(551.25 * eps, (double) 1 / (double) 9);
        h = numdiff::make_xph_representable(x[pos], h);

        auto x_h = new double[input_size];
        auto yh = new double[output_size];
        auto ymh = new double[output_size];
        auto y1 = new double[output_size];
        auto y2 = new double[output_size];
        auto y3 = new double[output_size];
        auto y4 = new double[output_size];
        auto tmp_y2 = new double[output_size];
        auto tmp_y3 = new double[output_size];
        auto tmp_y4 = new double[output_size];
        auto tmp1 = new double[output_size];
        auto tmp2 = new double[output_size];


        for (int i = 0; i < input_size; i++) x_h[i] = x[i];

        x_h[pos] = x[pos] + h;
        f(x_h, yh);

        x_h[pos] = x[pos] - h;
        f(x_h, ymh);

        x_h[pos] = x[pos] - 2 * h;
        f(x_h, y2);
        x_h[pos] = x[pos] + 2 * h;
        f(x_h, tmp_y2);

        x_h[pos] = x[pos] + 3 * h;
        f(x_h, y3);
        x_h[pos] = x[pos] - 3 * h;
        f(x_h, tmp_y3);

        x_h[pos] = x[pos] - 4 * h;
        f(x_h, y4);
        x_h[pos] = x[pos] + 4 * h;
        f(x_h, tmp_y4);

        for (int i = 0; i < output_size; i++) {
            y1[i] = yh[i] - ymh[i];
            y2[i] -= tmp_y2[i];
            y3[i] -= tmp_y3[i];
            y4[i] -= tmp_y4[i];
            tmp1[i] = 3 * y4[i] / 8 + 4 * y3[i];
            tmp2[i] = 21 * y2[i] + 84 * y1[i];
            y[i] = (tmp1[i] + tmp2[i]) / (105 * h);
        }

        delete[] x_h;
        delete[] yh;
        delete[] ymh;
        delete[] y1;
        delete[] y2;
        delete[] y3;
        delete[] y4;
        delete[] tmp_y2;
        delete[] tmp_y3;
        delete[] tmp_y4;
        delete[] tmp1;
        delete[] tmp2;

    }

}  // namespaces
#endif
