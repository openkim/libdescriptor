#include "Bispectrum.hpp"

#include <cmath>
#include <complex>
#include <limits>
#include <fstream>
//#include <numeric>

#ifdef MY_PI
#undef MY_PI
#endif

#define MY_PI 3.1415926535897932

#ifdef DIM
#undef DIM
#endif

/*! Problem dimensionality */
#define DIM 3

#ifdef LOG_ERROR
#undef LOG_ERROR
#endif

Bispectrum::Bispectrum(double const rfac0_in,
                       int const twojmax_in,
                       int const diagonalstyle_in,
                       int const use_shared_arrays_in,
                       double const rmin0_in,
                       int const switch_flag_in,
                       int const bzero_flag_in) :
        nmax(15),// TODO Fixed nmax, modify it after enzyme bugfix
        twojmax(twojmax_in),
        diagonalstyle(diagonalstyle_in),
        rmin0(rmin0_in),
        rfac0(rfac0_in),
        use_shared_arrays(use_shared_arrays_in),
        switch_flag(switch_flag_in),
        wself(1.0),
        bzero_flag(bzero_flag_in) {
    ncoeff = compute_ncoeff();

    create_twojmax_arrays();

    if (bzero_flag) {
        double const www = wself * wself * wself;
        for (int j = 1; j <= twojmax + 1; ++j) { bzero[j] = www * j; }
    }

    // 1D
    bvec.resize(ncoeff, static_cast<double>(0));

    // 2D
    dbvec.resize(ncoeff, 3);

    build_indexlist();

    init();

    grow_rij(nmax);
}

Bispectrum::~Bispectrum() {}

void Bispectrum::build_indexlist() {
    switch (diagonalstyle) {
        case 0: {
            int idxj_count = 0;
            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j2 = 0; j2 <= j1; ++j2) {
                    for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2);
                         j += 2) {
                        idxj_count++;
                    }
                }
            }

            // indexList can be changed here
            idxj.resize(idxj_count);

            idxj_max = idxj_count;

            idxj_count = 0;
            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j2 = 0; j2 <= j1; ++j2) {
                    for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2);
                         j += 2) {
                        idxj[idxj_count].j1 = j1;
                        idxj[idxj_count].j2 = j2;
                        idxj[idxj_count].j = j;
                        idxj_count++;
                    }
                }
            }
            return;
        }
        case 1: {
            int idxj_count = 0;
            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j = 0; j <= std::min(twojmax, 2 * j1); j += 2) {
                    idxj_count++;
                }
            }

            // indexList can be changed here
            idxj.resize(idxj_count);

            idxj_max = idxj_count;

            idxj_count = 0;
            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j = 0; j <= std::min(twojmax, 2 * j1); j += 2) {
                    idxj[idxj_count].j1 = j1;
                    idxj[idxj_count].j2 = j1;
                    idxj[idxj_count].j = j;
                    idxj_count++;
                }
            }
            return;
        }
        case 2: {
            int idxj_count = 0;
            for (int j1 = 0; j1 <= twojmax; ++j1) { idxj_count++; }

            // indexList can be changed here
            idxj.resize(idxj_count);
            idxj_max = idxj_count;
            idxj_count = 0;

            for (int j1 = 0; j1 <= twojmax; ++j1) {
                idxj[idxj_count].j1 = j1;
                idxj[idxj_count].j2 = j1;
                idxj[idxj_count].j = j1;
                idxj_count++;
            }
            return;
        }
        case 3: {
            int idxj_count = 0;
            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j2 = 0; j2 <= j1; ++j2) {
                    for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2);
                         j += 2) {
                        if (j >= j1) { idxj_count++; }
                    }
                }
            }

            // indexList can be changed here
            idxj.resize(idxj_count);
            idxj_max = idxj_count;
            idxj_count = 0;

            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j2 = 0; j2 <= j1; ++j2) {
                    for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2);
                         j += 2) {
                        if (j >= j1) {
                            idxj[idxj_count].j1 = j1;
                            idxj[idxj_count].j2 = j2;
                            idxj[idxj_count].j = j;
                            idxj_count++;
                        }
                    }
                }
            }
            return;
        }
        default:
            std::cerr << "The input style index = " + std::to_string(diagonalstyle)
                         + " is not a valid index!";
            std::abort();
    }
}

void Bispectrum::init() {
    init_clebsch_gordan();

    init_rootpqarray();
}

double Bispectrum::memory_usage() {
    double bytes;
    int const jdim = twojmax + 1;
    bytes = jdim * jdim * jdim * jdim * jdim * sizeof(double);
    bytes += 2 * jdim * jdim * jdim * sizeof(std::complex<double>);
    bytes += 2 * jdim * jdim * jdim * sizeof(double);
    bytes += jdim * jdim * jdim * 3 * sizeof(std::complex<double>);
    bytes += jdim * jdim * jdim * 3 * sizeof(double);
    bytes += ncoeff * sizeof(double);
    bytes += jdim * jdim * jdim * jdim * jdim * sizeof(std::complex<double>);
    return bytes;
}

//void Bispectrum::grow_rij(int const newnmax) {
//    if (newnmax > nmax) {
//        nmax = newnmax;
//
//        if (!use_shared_arrays) {
//                rij.resize(nmax * 3, 0.0);
//                inside.resize(nmax, 0);
//                wj.resize(nmax, 0.0);
//                rcutij.resize(nmax, 0.0);
//        }
//    }
//}

void Bispectrum::grow_rij(int const newnmax) {
    if (newnmax <= nmax) { return; }

    nmax = newnmax;

    if (!use_shared_arrays) {
        // 2D
        rij.resize(nmax, 3);

        // 1D
        inside.resize(nmax, 0);
        wj.resize(nmax, static_cast<double>(0));
        rcutij.resize(nmax, static_cast<double>(0));
    }
}

void Bispectrum::compute(int const index,
                         int const n_atoms,
                         int *const species,
                         int *const neigh_list,
                         int const number_of_neigh,
                         double *const coords,
                         double *const desc) {
    // prepare data
    Array2DView<double> coordinates(n_atoms, DIM, coords);

    int const *const ilist = neigh_list;
    int const iSpecies = species[index];

    // grow_rij removed for enzyme
    // insure rij, inside, wj, and rcutij are of size jnum

    // rij[, 3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    int ninside = 0;

    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < number_of_neigh; ++jj) {
        // adjust index of particle neighbor
        int const j = ilist[jj];
        int const jSpecies = species[j];

        // rij vec and
        double rvec[DIM];
        for (int dim = 0; dim < DIM; ++dim) {
            rvec[dim] = coordinates(j, dim) - coordinates(index, dim);
        }

        double const rsq
                = rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2];
        double const rmag = std::sqrt(rsq);

        // if (rmag < rcuts[iSpecies * n_species + jSpecies] && rmag > 1e-10) {
        if (rmag < rcuts(iSpecies, jSpecies) && rmag > 1e-10) {
            // rij[ninside * 3 + 0] = rvec[0];
            // rij[ninside * 3 + 1] = rvec[1];
            // rij[ninside * 3 + 2] = rvec[2];
            rij(ninside, 0) = rvec[0];
            rij(ninside, 1) = rvec[1];
            rij(ninside, 2) = rvec[2];

            inside[ninside] = j;
            wj[ninside] = wjelem[jSpecies];
            // rcutij[ninside] = rcuts[iSpecies * n_species  + jSpecies];
             rcutij[ninside] = rcuts(iSpecies, jSpecies);

            ninside++;
        }
    }

    // compute Ui, Zi, and Bi for atom I
    compute_ui(ninside);
    compute_zi();
    compute_bi();
    copy_bi2bvec();

    for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        desc[index * ncoeff + icoeff] = bvec[icoeff];
    }
}

void Bispectrum::set_cutoff(const char *name,
                            std::size_t const Nspecies,
                            double const *rcuts_in) {
    // store number of species and cutoff values

    for (int i = 0; i < Nspecies * Nspecies; i++) {
        rcuts.push_back(*(rcuts_in + i));
    }
}

void Bispectrum::set_weight(int const Nspecies, double const *weight_in) {
    wjelem.resize(Nspecies);
    std::copy(weight_in, weight_in + Nspecies, wjelem.data());
}

void Bispectrum::compute_ui(int const jnum) {
    zero_uarraytot();
    addself_uarraytot(wself);

    double x, y, z;
    double r, rsq;
    double z0;
    double theta0;

    for (int j = 0; j < jnum; ++j) {
        // x = rij[j * 3 + 0];
        // y = rij[j * 3 + 1];
        // z = rij[j * 3 + 2];

        x = rij(j, 0);
        y = rij(j, 1);
        z = rij(j, 2);

        rsq = x * x + y * y + z * z;
        r = std::sqrt(rsq);

        // TODO this is not in agreement with the paper, maybe change it
        theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij[j] - rmin0);
        // theta0 = (r - rmin0) * rscale0;

        z0 = r / std::tan(theta0);
        compute_uarray(x, y, z, z0, r);
        add_uarraytot(r, wj[j], rcutij[j]);
    }
}

void Bispectrum::compute_zi() {
    for (int j1 = 0; j1 <= twojmax; ++j1) {
        for (int j2 = 0; j2 <= j1; ++j2) {
            for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
                for (int mb = 0; 2 * mb <= j; ++mb) {
                    for (int ma = 0; ma <= j; ++ma) {
                        zarray_r(j1, j2, j, ma, mb) = 0.0;
                        zarray_i(j1, j2, j, ma, mb) = 0.0;

                        for (int ma1 = std::max(0, (2 * ma - j - j2 + j1) / 2);
                             ma1 <= std::min(j1, (2 * ma - j + j2 + j1) / 2);
                             ma1++) {
                            double sumb1_r = 0.0;
                            double sumb1_i = 0.0;

                            int const ma2 = (2 * ma - j - (2 * ma1 - j1) + j2) / 2;

                            for (int mb1 = std::max(0, (2 * mb - j - j2 + j1) / 2);
                                 mb1 <= std::min(j1, (2 * mb - j + j2 + j1) / 2);
                                 mb1++) {
                                int const mb2 = (2 * mb - j - (2 * mb1 - j1) + j2) / 2;

                                sumb1_r
                                        += cgarray(j1, j2, j, mb1, mb2)
                                           * (uarraytot_r(j1, ma1, mb1) * uarraytot_r(j2, ma2, mb2)
                                              - uarraytot_i(j1, ma1, mb1)
                                                * uarraytot_i(j2, ma2, mb2));
                                sumb1_i
                                        += cgarray(j1, j2, j, mb1, mb2)
                                           * (uarraytot_r(j1, ma1, mb1) * uarraytot_i(j2, ma2, mb2)
                                              + uarraytot_i(j1, ma1, mb1)
                                                * uarraytot_r(j2, ma2, mb2));
                            }

                            zarray_r(j1, j2, j, ma, mb)
                                    += sumb1_r * cgarray(j1, j2, j, ma1, ma2);
                            zarray_i(j1, j2, j, ma, mb)
                                    += sumb1_i * cgarray(j1, j2, j, ma1, ma2);
                        }
                    }
                }
            }
        }
    }
}

void Bispectrum::compute_bi() {
    for (int j1 = 0; j1 <= twojmax; ++j1) {
        for (int j2 = 0; j2 <= j1; ++j2) {
            for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2); j += 2) {
                barray(j1, j2, j) = 0.0;

                for (int mb = 0; 2 * mb < j; ++mb) {
                    for (int ma = 0; ma <= j; ++ma) {
                        barray(j1, j2, j)
                                += uarraytot_r(j, ma, mb) * zarray_r(j1, j2, j, ma, mb)
                                   + uarraytot_i(j, ma, mb) * zarray_i(j1, j2, j, ma, mb);
                    }
                }

                // For j even, special treatment for middle column
                if (j % 2 == 0) {
                    int const mb = j / 2;
                    for (int ma = 0; ma < mb; ++ma) {
                        barray(j1, j2, j)
                                += uarraytot_r(j, ma, mb) * zarray_r(j1, j2, j, ma, mb)
                                   + uarraytot_i(j, ma, mb) * zarray_i(j1, j2, j, ma, mb);
                    }

                    // ma = mb
                    barray(j1, j2, j)
                            += (uarraytot_r(j, mb, mb) * zarray_r(j1, j2, j, mb, mb)
                                + uarraytot_i(j, mb, mb) * zarray_i(j1, j2, j, mb, mb))
                               * 0.5;
                }

                barray(j1, j2, j) *= 2.0;

                if (bzero_flag) { barray(j1, j2, j) -= bzero[j]; }
            }
        }
    }
}

void Bispectrum::copy_bi2bvec() {
    switch (diagonalstyle) {
        case (0): {
            for (int j1 = 0, ncount = 0; j1 <= twojmax; ++j1) {
                for (int j2 = 0; j2 <= j1; ++j2) {
                    for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2);
                         j += 2) {
                        bvec[ncount++] = barray(j1, j2, j);
                    }
                }
            }
            return;
        }
        case (1): {
            for (int j1 = 0, ncount = 0; j1 <= twojmax; ++j1) {
                for (int j = 0; j <= std::min(twojmax, 2 * j1); j += 2) {
                    bvec[ncount++] = barray(j1, j1, j);
                }
            }
            return;
        }
        case (2): {
            for (int j1 = 0, ncount = 0; j1 <= twojmax; ++j1) {
                bvec[ncount++] = barray(j1, j1, j1);
            }
            return;
        }
        case (3): {
            for (int j1 = 0, ncount = 0; j1 <= twojmax; ++j1) {
                for (int j2 = 0; j2 <= j1; ++j2) {
                    for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2);
                         j += 2) {
                        if (j >= j1) { bvec[ncount++] = barray(j1, j2, j); }
                    }
                }
            }
            return;
        }
    }
}

void Bispectrum::zero_uarraytot() {
    for (int j = 0; j <= twojmax; ++j) {
        for (int ma = 0; ma <= j; ++ma) {
            for (int mb = 0; mb <= j; ++mb) {
                uarraytot_r(j, ma, mb) = 0.0;
                uarraytot_i(j, ma, mb) = 0.0;
            }
        }
    }
}

void Bispectrum::addself_uarraytot(double const wself_in) {
    for (int j = 0; j <= twojmax; ++j) {
        for (int ma = 0; ma <= j; ++ma) {
            uarraytot_r(j, ma, ma) = wself_in;
            uarraytot_i(j, ma, ma) = 0.0;
        }
    }
}

void Bispectrum::compute_uarray(double const x,
                                double const y,
                                double const z,
                                double const z0,
                                double const r) {
    // compute Cayley-Klein parameters for unit quaternion
    double const r0inv = 1.0 / std::sqrt(r * r + z0 * z0);

    double const a_r = r0inv * z0;
    double const a_i = -r0inv * z;
    double const b_r = r0inv * y;
    double const b_i = -r0inv * x;

    double rootpq;

    // VMK Section 4.8.2
    uarray_r(0, 0, 0) = 1.0;
    uarray_i(0, 0, 0) = 0.0;

    for (int j = 1; j <= twojmax; ++j) {
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2 * mb <= j; ++mb) {
            uarray_r(j, 0, mb) = 0.0;
            uarray_i(j, 0, mb) = 0.0;

            for (int ma = 0; ma < j; ++ma) {
                rootpq = rootpqarray[(j - ma) * twojmax + (j - mb)];
                uarray_r(j, ma, mb) += rootpq
                                       * (a_r * uarray_r(j - 1, ma, mb)
                                          + a_i * uarray_i(j - 1, ma, mb));
                uarray_i(j, ma, mb) += rootpq
                                       * (a_r * uarray_i(j - 1, ma, mb)
                                          - a_i * uarray_r(j - 1, ma, mb));

                rootpq = rootpqarray[(ma + 1) * twojmax + (j - mb)];
                uarray_r(j, ma + 1, mb)
                        = -rootpq
                          * (b_r * uarray_r(j - 1, ma, mb) + b_i * uarray_i(j - 1, ma, mb));
                uarray_i(j, ma + 1, mb)
                        = -rootpq
                          * (b_r * uarray_i(j - 1, ma, mb) - b_i * uarray_r(j - 1, ma, mb));
            }
        }

        // copy left side to right side with inversion symmetry VMK 4.4(2)
        // u[ma-j, mb-j] = (-1)^(ma-mb)*Conj([u[ma, mb])

        int mbpar = -1;
        for (int mb = 0; 2 * mb <= j; ++mb) {
            mbpar = -mbpar;

            int mapar = -mbpar;
            for (int ma = 0; ma <= j; ++ma) {
                mapar = -mapar;
                if (mapar == 1) {
                    uarray_r(j, j - ma, j - mb) = uarray_r(j, ma, mb);
                    uarray_i(j, j - ma, j - mb) = -uarray_i(j, ma, mb);
                } else {
                    uarray_r(j, j - ma, j - mb) = -uarray_r(j, ma, mb);
                    uarray_i(j, j - ma, j - mb) = uarray_i(j, ma, mb);
                }
            }
        }
    }
}

void Bispectrum::add_uarraytot(double const r,
                               double const wj_in,
                               double const rcut_in) {
    double sfac = compute_sfac(r, rcut_in);
    sfac *= wj_in;

    for (int j = 0; j <= twojmax; ++j) {
        for (int ma = 0; ma <= j; ++ma) {
            for (int mb = 0; mb <= j; ++mb) {
                uarraytot_r(j, ma, mb) += sfac * uarray_r(j, ma, mb);
                uarraytot_i(j, ma, mb) += sfac * uarray_i(j, ma, mb);
            }
        }
    }
}

double Bispectrum::compute_sfac(double const r, double const rcut_in) {
    switch (switch_flag) {
        case (0):
            return 1.0;
        case (1):
            return (r <= rmin0) ? 1.0
                                : (r > rcut_in)
                                  ? 0.0
                                  : 0.5
                                    * (std::cos((r - rmin0) * MY_PI / (rcut_in - rmin0))
                                       + 1.0);
        default:
            return 0.0;
    }
}

double Bispectrum::compute_dsfac(double const r, double const rcut_in) {
    switch (switch_flag) {
        case (1):
            if (r <= rmin0 || r > rcut_in) { return 0.0; }
            else {
                double const rcutfac = MY_PI / (rcut_in - rmin0);
                return -0.5 * std::sin((r - rmin0) * rcutfac) * rcutfac;
            }
        default:
            return 0.0;
    }
}

void Bispectrum::create_twojmax_arrays() {
    int const jdim = twojmax + 1;

    cgarray.resize(jdim, jdim, jdim, jdim, jdim);
    rootpqarray.resize((jdim + 1) * (jdim + 1));
    barray.resize(jdim, jdim, jdim);
    dbarray.resize(jdim, jdim, jdim, 3);
    duarray_r.resize(jdim, jdim, jdim, 3);
    duarray_i.resize(jdim, jdim, jdim, 3);
    uarray_r.resize(jdim, jdim, jdim);
    uarray_i.resize(jdim, jdim, jdim);

    if (bzero_flag) { bzero.resize(jdim, static_cast<double>(0)); }

    if (!use_shared_arrays) {
        uarraytot_r.resize(jdim, jdim, jdim);
        uarraytot_i.resize(jdim, jdim, jdim);

        zarray_r.resize(jdim, jdim, jdim, jdim, jdim);
        zarray_i.resize(jdim, jdim, jdim, jdim, jdim);
    }
}

inline double Bispectrum::factorial(int const n) {
    if (n < 0 || n > nmaxfactorial) {
        std::cerr << "The input n = " + std::to_string(n)
                     + " is not valid for a factorial!";
        std::abort();
    }
    return nfac_table[n];
}

double const Bispectrum::nfac_table[] = {
        1,
        1,
        2,
        6,
        24,
        120,
        720,
        5040,
        40320,
        362880,
        3628800,
        39916800,
        479001600,
        6227020800,
        87178291200,
        1307674368000,
        20922789888000,
        355687428096000,
        6.402373705728e+15,
        1.21645100408832e+17,
        2.43290200817664e+18,
        5.10909421717094e+19,
        1.12400072777761e+21,
        2.5852016738885e+22,
        6.20448401733239e+23,
        1.5511210043331e+25,
        4.03291461126606e+26,
        1.08888694504184e+28,
        3.04888344611714e+29,
        8.8417619937397e+30,
        2.65252859812191e+32,
        8.22283865417792e+33,
        2.63130836933694e+35,
        8.68331761881189e+36,
        2.95232799039604e+38,
        1.03331479663861e+40,
        3.71993326789901e+41,
        1.37637530912263e+43,
        5.23022617466601e+44,
        2.03978820811974e+46,
        8.15915283247898e+47,
        3.34525266131638e+49,
        1.40500611775288e+51,
        6.04152630633738e+52,
        2.65827157478845e+54,
        1.1962222086548e+56,
        5.50262215981209e+57,
        2.58623241511168e+59,
        1.24139155925361e+61,
        6.08281864034268e+62,
        3.04140932017134e+64,
        1.55111875328738e+66,
        8.06581751709439e+67,
        4.27488328406003e+69,
        2.30843697339241e+71,
        1.26964033536583e+73,
        7.10998587804863e+74,
        4.05269195048772e+76,
        2.35056133128288e+78,
        1.3868311854569e+80,
        8.32098711274139e+81,
        5.07580213877225e+83,
        3.14699732603879e+85,
        1.98260831540444e+87,
        1.26886932185884e+89,
        8.24765059208247e+90,
        5.44344939077443e+92,
        3.64711109181887e+94,
        2.48003554243683e+96,
        1.71122452428141e+98,
        1.19785716699699e+100,
        8.50478588567862e+101,
        6.12344583768861e+103,
        4.47011546151268e+105,
        3.30788544151939e+107,
        2.48091408113954e+109,
        1.88549470166605e+111,
        1.45183092028286e+113,
        1.13242811782063e+115,
        8.94618213078297e+116,
        7.15694570462638e+118,
        5.79712602074737e+120,
        4.75364333701284e+122,
        3.94552396972066e+124,
        3.31424013456535e+126,
        2.81710411438055e+128,
        2.42270953836727e+130,
        2.10775729837953e+132,
        1.85482642257398e+134,
        1.65079551609085e+136,
        1.48571596448176e+138,
        1.3520015276784e+140,
        1.24384140546413e+142,
        1.15677250708164e+144,
        1.08736615665674e+146,
        1.03299784882391e+148,
        9.91677934870949e+149,
        9.61927596824821e+151,
        9.42689044888324e+153,
        9.33262154439441e+155,
        9.33262154439441e+157,
        9.42594775983835e+159,
        9.61446671503512e+161,
        9.90290071648618e+163,
        1.02990167451456e+166,
        1.08139675824029e+168,
        1.14628056373471e+170,
        1.22652020319614e+172,
        1.32464181945183e+174,
        1.44385958320249e+176,
        1.58824554152274e+178,
        1.76295255109024e+180,
        1.97450685722107e+182,
        2.23119274865981e+184,
        2.54355973347219e+186,
        2.92509369349301e+188,
        3.3931086844519e+190,
        3.96993716080872e+192,
        4.68452584975429e+194,
        5.5745857612076e+196,
        6.68950291344912e+198,
        8.09429852527344e+200,
        9.8750442008336e+202,
        1.21463043670253e+205,
        1.50614174151114e+207,
        1.88267717688893e+209,
        2.37217324288005e+211,
        3.01266001845766e+213,
        3.8562048236258e+215,
        4.97450422247729e+217,
        6.46685548922047e+219,
        8.47158069087882e+221,
        1.118248651196e+224,
        1.48727070609069e+226,
        1.99294274616152e+228,
        2.69047270731805e+230,
        3.65904288195255e+232,
        5.01288874827499e+234,
        6.91778647261949e+236,
        9.61572319694109e+238,
        1.34620124757175e+241,
        1.89814375907617e+243,
        2.69536413788816e+245,
        3.85437071718007e+247,
        5.5502938327393e+249,
        8.04792605747199e+251,
        1.17499720439091e+254,
        1.72724589045464e+256,
        2.55632391787286e+258,
        3.80892263763057e+260,
        5.71338395644585e+262,
        8.62720977423323e+264,
        1.31133588568345e+267,
        2.00634390509568e+269,
        3.08976961384735e+271,
        4.78914290146339e+273,
        7.47106292628289e+275,
        1.17295687942641e+278,
        1.85327186949373e+280,
        2.94670227249504e+282,
        4.71472363599206e+284,
        7.59070505394721e+286,
        1.22969421873945e+289,
        2.0044015765453e+291,
        3.28721858553429e+293,
        5.42391066613159e+295,
        9.00369170577843e+297,
        1.503616514865e+300,  // nmaxfactorial = 167
};

inline double Bispectrum::deltacg(int const j1, int const j2, int const j) {
    double const sfaccg = factorial((j1 + j2 + j) / 2 + 1);
    return std::sqrt(factorial((j1 + j2 - j) / 2) * factorial((j1 - j2 + j) / 2)
                     * factorial((-j1 + j2 + j) / 2) / sfaccg);
}

void Bispectrum::init_clebsch_gordan() {
    for (int j1 = 0; j1 <= twojmax; ++j1) {
        for (int j2 = 0; j2 <= twojmax; ++j2) {
            for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2); j += 2) {
                for (int m1 = 0; m1 <= j1; m1 += 1) {
                    int const aa2 = 2 * m1 - j1;

                    for (int m2 = 0; m2 <= j2; m2 += 1) {
                        // -c <= cc <= c
                        int const bb2 = 2 * m2 - j2;
                        int const m = (aa2 + bb2 + j) / 2;
                        if (m < 0 || m > j) { continue; }

                        double sum(0.0);
                        for (int z = std::max(
                                0, std::max(-(j - j2 + aa2) / 2, -(j - j1 - bb2) / 2));
                             z <= std::min((j1 + j2 - j) / 2,
                                           std::min((j1 - aa2) / 2, (j2 + bb2) / 2));
                             z++) {
                            int const ifac = z % 2 ? -1 : 1;
                            sum += ifac
                                   / (factorial(z) * factorial((j1 + j2 - j) / 2 - z)
                                      * factorial((j1 - aa2) / 2 - z)
                                      * factorial((j2 + bb2) / 2 - z)
                                      * factorial((j - j2 + aa2) / 2 + z)
                                      * factorial((j - j1 - bb2) / 2 + z));
                        }

                        int const cc2 = 2 * m - j;

                        double dcg = deltacg(j1, j2, j);

                        double sfaccg = std::sqrt(
                                factorial((j1 + aa2) / 2) * factorial((j1 - aa2) / 2)
                                * factorial((j2 + bb2) / 2) * factorial((j2 - bb2) / 2)
                                * factorial((j + cc2) / 2) * factorial((j - cc2) / 2)
                                * (j + 1));

                        cgarray(j1, j2, j, m1, m2) = sum * dcg * sfaccg;
                    }
                }
            }
        }
    }
}

void Bispectrum::init_rootpqarray() {
    for (int p = 1; p <= twojmax; p++) {
        for (int q = 1; q <= twojmax; q++) {
            rootpqarray[p * twojmax + q] = std::sqrt(static_cast<double>(p) / q);
        }
    }
}

inline void Bispectrum::jtostr(char *str_out, int const j) {
    if (j % 2 == 0) { sprintf(str_out, "%d", j / 2); }
    else {
        sprintf(str_out, "%d/2", j);
    }
}

inline void Bispectrum::mtostr(char *str_out, int const j, int const m) {
    if (j % 2 == 0) { sprintf(str_out, "%d", m - j / 2); }
    else {
        sprintf(str_out, "%d/2", 2 * m - j);
    }
}

void Bispectrum::print_clebsch_gordan(FILE *file) {
    char stra[20];
    char strb[20];
    char strc[20];
    char straa[20];
    char strbb[20];
    char strcc[20];

    int m, aa2, bb2;

    fprintf(file, "a, aa, b, bb, c, cc, c(a,aa,b,bb,c,cc) \n");

    for (int j1 = 0; j1 <= twojmax; ++j1) {
        jtostr(stra, j1);

        for (int j2 = 0; j2 <= twojmax; ++j2) {
            jtostr(strb, j2);

            for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2); j += 2) {
                jtostr(strc, j);

                for (int m1 = 0; m1 <= j1; ++m1) {
                    mtostr(straa, j1, m1);

                    aa2 = 2 * m1 - j1;

                    for (int m2 = 0; m2 <= j2; ++m2) {
                        bb2 = 2 * m2 - j2;
                        m = (aa2 + bb2 + j) / 2;

                        if (m < 0 || m > j) { continue; }

                        mtostr(strbb, j2, m2);

                        mtostr(strcc, j, m);

                        fprintf(file,
                                "%s\t%s\t%s\t%s\t%s\t%s\t%g\n",
                                stra,
                                straa,
                                strb,
                                strbb,
                                strc,
                                strcc,
                                cgarray(j1, j2, j, m1, m2));
                    }
                }
            }
        }
    }
}

int Bispectrum::compute_ncoeff() {
    switch (diagonalstyle) {
        case (0): {
            int ncount(0);
            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j2 = 0; j2 <= j1; ++j2) {
                    for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2);
                         j += 2) {
                        ncount++;
                    }
                }
            }
            return ncount;
        }
        case (1): {
            int ncount(0);
            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j = 0; j <= std::min(twojmax, 2 * j1); j += 2) { ncount++; }
            }
            return ncount;
        }
        case (2): {
            int ncount(0);
            for (int j1 = 0; j1 <= twojmax; ++j1) { ncount++; }
            return ncount;
        }
        case (3): {
            int ncount(0);
            for (int j1 = 0; j1 <= twojmax; ++j1) {
                for (int j2 = 0; j2 <= j1; ++j2) {
                    for (int j = std::abs(j1 - j2); j <= std::min(twojmax, j1 + j2);
                         j += 2) {
                        if (j >= j1) { ncount++; }
                    }
                }
            }
            return ncount;
        }
        default:
            std::cerr << "The input style index = " + std::to_string(diagonalstyle)
                         + " is not a valid index!!";
            std::abort();
    }
}

Bispectrum::Bispectrum(std::string &file_name) {
    initFromFile(file_name);
}

void Bispectrum::initFromFile(std::string &file_name) {
    std::fstream file_ptr(file_name);
    std::string placeholder_string;

    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    n_species = std::stoi(placeholder_string);
    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');

    auto cutoff_matrix = new double[n_species * n_species];
    for (int i = 0; i < n_species; i++) {
        for (int j = 0; j < n_species; j++) {
            auto pos = placeholder_string.find(' ');
            *(cutoff_matrix + n_species * i + j) = std::stod(placeholder_string.substr(0, pos));
            if (pos != std::string::npos) placeholder_string.erase(0, pos + 1);
        }
        std::getline(file_ptr, placeholder_string);
    }

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    std::string cutoff_function = placeholder_string;

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    width = std::stoi(placeholder_string);

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    int jmax = stoi(placeholder_string);

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    rfac0 = std::stod(placeholder_string);

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    diagonalstyle = std::stoi(placeholder_string);

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    rmin0 = std::stoi(placeholder_string);

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    switch_flag = std::stoi(placeholder_string);

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    bzero_flag = std::stoi(placeholder_string);

    auto weights = new double[n_species];
    for (int i = 0; i < n_species; i++){
        *(weights + i) = 1.0;
    }

    set_weight(n_species, weights);
    set_cutoff(cutoff_function.c_str(), n_species, cutoff_matrix);

    use_shared_arrays = 0;
    nmax = 25; // This value is maximum number of neighbors possible,
    // Ideally it should be determined at runtime, but enzyme has issues with it
    // so fixing the upperbound. Would be made an hyper-parameter soon, till
    // Enzyme people fixes it. TODO
    grow_rij(nmax);

    twojmax = 2 * jmax;
    wself = 1.0;

    ncoeff = compute_ncoeff();

    create_twojmax_arrays();

    if (bzero_flag) {
        double const www = wself * wself * wself;
        for (int j = 1; j <= twojmax + 1; ++j) { bzero[j] = www * j; }
    }

    // 1D
    bvec.resize(ncoeff, static_cast<double>(0));

    // 2D
    dbvec.resize(ncoeff, 3);

    build_indexlist();

    init();

    delete[] cutoff_matrix;
    delete[] weights;
}

void Bispectrum::clone_empty(DescriptorKind *ds) {}

#undef MY_PI
#undef DIM
