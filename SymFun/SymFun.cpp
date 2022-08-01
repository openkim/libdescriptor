#include "SymFun.hpp"
#include <vector>
#include <cmath>

#ifdef DIM
#undef DIM
#endif
#define DIM 3

#ifdef MY_PI
#undef MY_PI
#endif
#define MY_PI 3.1415926535897932

int enzyme_dup;
int enzyme_out;
int enzyme_const;

inline double cut_cos(double const r, double const rcut) {
    return (r < rcut) ? 0.5 * (std::cos(MY_PI * r / rcut) + 1.0) : 0.0;
}


void sym_g5(double const zeta, double const lambda, double const eta, double const *r,
            double const *rcut, double &phi) {
    double const rij = r[0];
    double const rik = r[1];
    double const rcutij = rcut[0];
    double const rcutik = rcut[1];

    if (rij > rcutij || rik > rcutik) { phi = 0.0; }
    else {
        double const rjk = r[2];
        double const rijsq = rij * rij;
        double const riksq = rik * rik;
        double const rjksq = rjk * rjk;

        // i is the apex atom
        double const cos_ijk = (rijsq + riksq - rjksq) / (2 * rij * rik);
        double const base = 1.0 + lambda * cos_ijk;

        // prevent numerical instability (when lambda=-1 and cos_ijk=1)
        double const costerm = (base <= 0) ? 0.0 : std::pow(base, zeta);
        double const eterm = std::exp(-eta * (rijsq + riksq));
        phi = std::pow(2, 1 - zeta) * costerm * eterm * cut_cos(rij, rcutij)
              * cut_cos(rik, rcutik);
    }
}

void sym_g4(double const zeta, double const lambda, double const eta, double const *r,
            double const *rcut, double &phi) {
    double const rij = r[0];
    double const rik = r[1];
    double const rjk = r[2];
    double const rcutij = rcut[0];
    double const rcutik = rcut[1];
    double const rcutjk = rcut[2];

    if (rij > rcutij || rik > rcutik || rjk > rcutjk) { phi = 0.0; }
    else {
        double const rijsq = rij * rij;
        double const riksq = rik * rik;
        double const rjksq = rjk * rjk;

        // i is the apex atom
        double const cos_ijk = (rijsq + riksq - rjksq) / (2 * rij * rik);
        double const base = 1 + lambda * cos_ijk;

        // prevent numerical instability (when lambda=-1 and cos_ijk=1)
        double const costerm = (base <= 0) ? 0.0 : std::pow(base, zeta);
        double const eterm = std::exp(-eta * (rijsq + riksq + rjksq));

        phi = std::pow(2, 1 - zeta) * costerm * eterm * cut_cos(rij, rcutij)
              * cut_cos(rik, rcutik) * cut_cos(rjk, rcutjk);
    }
}

void sym_g3(double const kappa, double const r, double const rcut, double &phi) {
    phi = std::cos(kappa * r) * cut_cos(r, rcut);
}

void sym_g2(double const eta, double const Rs, double const r, double const rcut, double &phi) {
    phi = std::exp(-eta * (r - Rs) * (r - Rs)) * cut_cos(r, rcut);
}

void sym_g1(double const r, double const rcut, double &phi) {
    phi = cut_cos(r, rcut);
}


void symmetry_function_atomic(int const i,
                              double const *coords,
                              int const *particleSpeciesCodes,
                              int const *neighlist,
                              int const numnei,
                              double *const desc,
                              SymmetryFunctionParams *SymParam) {
    // prepare data
    VectorOfSizeDIM *coordinates = (VectorOfSizeDIM *) coords;
    int const iSpecies = particleSpeciesCodes[i];
    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < numnei; ++jj) {
        // adjust index of particle neighbor
        int const j = neighlist[jj];
        int const jSpecies = particleSpeciesCodes[j];
        // cutoff between ij
        double rcutij = SymParam->rcut_2D_(iSpecies, jSpecies);
        // Compute rij
        double rij[DIM];
        for (int dim = 0; dim < DIM; ++dim) {
            rij[dim] = coordinates[j][dim] - coordinates[i][dim];
        }

        double const rijsq = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
        double const rijmag = std::sqrt(rijsq);

        // if particles i and j not interact
        if (rijmag > rcutij) { continue; }

        // Loop over descriptors
        // two-body descriptors
        for (std::size_t p = 0; p < SymParam->name_.size(); ++p) {
            if (SymParam->name_[p] != 1 && SymParam->name_[p] != 2 && SymParam->name_[p] != 3) {
                continue;
            }

            int idx = SymParam->starting_index_[p];
            // Loop over same descriptor but different parameter set
            for (int q = 0; q < SymParam->num_param_sets_[p]; ++q) {
                double gc = 0.0;

                if (SymParam->name_[p] == 1) {
                    sym_g1(rijmag, rcutij, gc);
                } else if (SymParam->name_[p] == 2) {
                    double eta = SymParam->params_[p](q, 0);
                    auto Rs = SymParam->params_[p](q, 1);
                    sym_g2(eta, Rs, rijmag, rcutij, gc);
                } else if (SymParam->name_[p] == 3) {
                    double kappa = SymParam->params_[p](q, 0);
                    sym_g3(kappa, rijmag, rcutij, gc);
                }

                desc[idx] += gc;
                ++idx;
            }
        }
        // three-body descriptors
        if (SymParam->has_three_body_ == 0) { continue; }

        // Loop over kk
        for (int kk = jj + 1; kk < numnei; ++kk) {
            // Adjust index of particle neighbor
            int const k = neighlist[kk];
            int const kSpecies = particleSpeciesCodes[k];

            // cutoff between ik and jk
            double const rcutik = SymParam->rcut_2D_[iSpecies][kSpecies];
            double const rcutjk = SymParam->rcut_2D_[jSpecies][kSpecies];

            // Compute rik, rjk and their squares
            double rik[DIM];
            double rjk[DIM];

            for (int dim = 0; dim < DIM; ++dim) {
                rik[dim] = coordinates[k][dim] - coordinates[i][dim];
                rjk[dim] = coordinates[k][dim] - coordinates[j][dim];
            }

            double const riksq = rik[0] * rik[0] + rik[1] * rik[1] + rik[2] * rik[2];
            double const rjksq = rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2];
            double const rikmag = std::sqrt(riksq);
            double const rjkmag = std::sqrt(rjksq);

            // Check whether three-dody not interacting
            if (rikmag > rcutik) { continue; }

            double const rvec[3] = {rijmag, rikmag, rjkmag};
            double const rcutvec[3] = {rcutij, rcutik, rcutjk};

            // Loop over descriptors
            // three-body descriptors
            for (size_t p = 0; p < SymParam->name_.size(); ++p) {
                if (SymParam->name_[p] != 4 && SymParam->name_[p] != 5) { continue; }
                int idx = SymParam->starting_index_[p];

                // Loop over same descriptor but different parameter set
                for (int q = 0; q < SymParam->num_param_sets_[p]; ++q) {
                    double gc = 0.0;

                    if (SymParam->name_[p] == 4) {
                        double zeta = SymParam->params_[p](q, 0);
                        double lambda = SymParam->params_[p](q, 1);
                        double eta = SymParam->params_[p](q, 2);

                        sym_g4(zeta, lambda, eta, rvec, rcutvec, gc);
                    } else if (SymParam->name_[p] == 5) {
                        double zeta = SymParam->params_[p](q, 0);
                        double lambda = SymParam->params_[p](q, 1);
                        double eta = SymParam->params_[p](q, 2);

                        sym_g5(zeta, lambda, eta, rvec, rcutvec, gc);
                    }

                    desc[idx] += gc;
                    ++idx;
                }
            }
        }
    }
}

void grad_symmetry_function_atomic(int const i,
                                   double const *coords,
                                   double const *d_coords,
                                   int const *particleSpeciesCodes,
                                   int const *neighlist,
                                   int const numnei,
                                   double *const desc,
                                   double const *d_grad_loss_zeta,
                                   SymmetryFunctionParams *SymParam) {
    __enzyme_autodiff(symmetry_function_atomic,
                      enzyme_const, i,
                      enzyme_dup, coords, d_coords,
                      enzyme_const, particleSpeciesCodes,
                      enzyme_const, neighlist,
                      enzyme_const, numnei,
                      enzyme_dup, desc, d_grad_loss_zeta,
                      enzyme_const, SymParam);
}