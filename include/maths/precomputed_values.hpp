#ifndef __PRECOMPUTED_VALUES_HPP__
#define __PRECOMPUTED_VALUES_HPP__

namespace precomputed_values {
        int const bessel_zeros_max_j = 100;
        int const bessel_zeros_max_roots = 100;
        extern double bessel_zeros[bessel_zeros_max_j + 1][bessel_zeros_max_roots];
}

#endif //__PRECOMPUTED_VALUES_HPP__