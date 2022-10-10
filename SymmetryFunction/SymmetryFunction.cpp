#include "SymmetryFunction.hpp"
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

inline void SymmetryFunctions::set_species(std::vector<std::string> &species) {
    species_.resize(species.size());
    std::copy(species.begin(), species.end(), species_.begin());
};

inline void SymmetryFunctions::get_species(std::vector<std::string> &species) {
    species.resize(species_.size());
    std::copy(species_.begin(), species_.end(), species.begin());
};

inline int SymmetryFunctions::get_num_species() { return species_.size(); }

void SymmetryFunctions::set_cutoff(char const *name,
                                   std::size_t const Nspecies,
                                   double const *rcut_2D) {
    (void) name;   // to avoid unused warning
    rcut_2D_.resize(Nspecies, Nspecies, rcut_2D);
}

inline double SymmetryFunctions::get_cutoff(int const iCode, int const jCode) {
    return rcut_2D_(iCode, jCode);
};

void SymmetryFunctions::add_descriptor(char const *name,
                                       double const *values,
                                       int const row,
                                       int const col) {
    if (strcmp(name, "g1") == 0) { name_.push_back(1); };
    if (strcmp(name, "g2") == 0) { name_.push_back(2); };
    if (strcmp(name, "g3") == 0) { name_.push_back(3); };
    if (strcmp(name, "g4") == 0) { name_.push_back(4); };
    if (strcmp(name, "g5") == 0) { name_.push_back(5); };

    Array2D<double> params(row, col, values);
    params_.push_back(std::move(params));

    auto sum = std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
    starting_index_.push_back(sum);

    num_param_sets_.push_back(row);
    num_params_.push_back(col);

    if (strcmp(name, "g4") == 0 || strcmp(name, "g5") == 0) {
        has_three_body_ = true;
    }
}

inline double cut_cos(double const r, double const rcut) {
    return (r < rcut) ? 0.5 * (std::cos(MY_PI * r / rcut) + 1.0) : 0.0;
}

int SymmetryFunctions::get_num_descriptors() {
    return std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
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

        // index is the apex atom
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

        // index is the apex atom
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

void SymmetryFunctions::compute(int const index,
                                int const n_atoms,
                                int *const species,
                                int *const neigh_list,
                                int const number_of_neigh,
                                double *const coords,
                                double *const desc) {
    // prepare data
    auto *coordinates = (VectorOfSizeDIM *) coords;
    int const iSpecies = species[index];
    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < number_of_neigh; ++jj) {
        // adjust index of particle neighbor
        int const j = neigh_list[jj];
        int const jSpecies = species[j];
        // cutoff between ij
        double rcutij = rcut_2D_(iSpecies, jSpecies);
        // Compute rij
        double rij[DIM];
        for (int dim = 0; dim < DIM; ++dim) {
            rij[dim] = coordinates[j][dim] - coordinates[index][dim];
        }

        double const rijsq = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
        double const rijmag = std::sqrt(rijsq);

        // if particles index and j not interact
        if (rijmag > rcutij) { continue; }

        // Loop over descriptors
        // two-body descriptors
        for (std::size_t p = 0; p < name_.size(); ++p) {
            if (name_[p] != 1 && name_[p] != 2 && name_[p] != 3) {
                continue;
            }

            int idx = starting_index_[p];
            // Loop over same descriptor but different parameter set
            for (int q = 0; q < num_param_sets_[p]; ++q) {
                double gc = 0.0;

                if (name_[p] == 1) {
                    sym_g1(rijmag, rcutij, gc);
                } else if (name_[p] == 2) {
                    double eta = params_[p](q, 0);
                    auto Rs = params_[p](q, 1);
                    sym_g2(eta, Rs, rijmag, rcutij, gc);
                } else if (name_[p] == 3) {
                    double kappa = params_[p](q, 0);
                    sym_g3(kappa, rijmag, rcutij, gc);
                }
                desc[idx] += gc;
                ++idx;
            }
        }
        // three-body descriptors
        if (has_three_body_ == 0) { continue; }

        // Loop over kk
        for (int kk = jj + 1; kk < number_of_neigh; ++kk) {
            // Adjust index of particle neighbor
            int const k = neigh_list[kk];
            int const kSpecies = species[k];

            // cutoff between ik and jk
            double const rcutik = rcut_2D_[iSpecies][kSpecies];
            double const rcutjk = rcut_2D_[jSpecies][kSpecies];

            // Compute rik, rjk and their squares
            double rik[DIM];
            double rjk[DIM];

            for (int dim = 0; dim < DIM; ++dim) {
                rik[dim] = coordinates[k][dim] - coordinates[index][dim];
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
            for (size_t p = 0; p < name_.size(); ++p) {
                if (name_[p] != 4 && name_[p] != 5) { continue; }
                int idx = starting_index_[p];

                // Loop over same descriptor but different parameter set
                for (int q = 0; q < num_param_sets_[p]; ++q) {
                    double gc = 0.0;

                    if (name_[p] == 4) {
                        double zeta = params_[p](q, 0);
                        double lambda = params_[p](q, 1);
                        double eta = params_[p](q, 2);

                        sym_g4(zeta, lambda, eta, rvec, rcutvec, gc);
                    } else if (name_[p] == 5) {
                        double zeta = params_[p](q, 0);
                        double lambda = params_[p](q, 1);
                        double eta = params_[p](q, 2);

                        sym_g5(zeta, lambda, eta, rvec, rcutvec, gc);
                    }

                    desc[idx] += gc;
                    ++idx;
                }
            }
        }
    }
}

SymmetryFunctions::SymmetryFunctions(std::string &file_name) {

};