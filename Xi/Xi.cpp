#include "Xi.hpp"
#include "get_parameters.hpp"
#include <cmath>
# include "maths/maths.hpp"

#define MAX_NEIGHBORS 100

Xi::Xi(int q, double cutoff, std::vector<std::string> &species, std::string &radial_basis) {
    this->cutoff = cutoff;
    this->radial_basis = radial_basis;
    this->q = q;
    this->species_ = species;
    this->width = get_width();
    this->ln_params.resize(this->width * 6);
    copy_params(this->width, this->ln_params.data()); // get the parameters from the file
    this->l_max_sq = (q + 1) * (q + 1);
    if (radial_basis != "bessel") {
        throw std::runtime_error("Radial basis " + radial_basis + " not implemented");
    }
//    this->create_clebsh_gordon_array();
}

void Xi::create_clebsh_gordon_array() {
}

int Xi::get_width() {
    return this->width;
}

void Xi::compute(int index,
                 int n_atoms,
                 int *species,
                 int *neighbor_lists,
                 int number_of_neighbors,
                 double *coordinates,
                 double *desc) {
    auto *coordinates_  = (VectorOfSize3 *) coordinates;
    std::array<double, 3> r_i = {coordinates_[index][0], coordinates_[index][1], coordinates_[index][2]};

    auto center_shifted_neighbor_coordinates = std::vector<double>(3 * MAX_NEIGHBORS);
    auto center_shifted_neighbor_coordinates_zj = std::vector<double>(3 * MAX_NEIGHBORS);

    //1. get center shifted coordinates for index atom. Neighors exclude the center atom
    auto i_coordinates = (VectorOfSize3 *) center_shifted_neighbor_coordinates.data();
    bool rigid_shift = false;
    for (int i = 0; i < number_of_neighbors; i++) {
        int neighbor_index = neighbor_lists[i];
        std::array<double, 3> r_j = {coordinates_[neighbor_index][0], coordinates_[neighbor_index][1],
                                     coordinates_[neighbor_index][2]};
        i_coordinates[i][0] = r_j[0] - r_i[0];
        i_coordinates[i][1] = r_j[1] - r_i[1];
        i_coordinates[i][2] = r_j[2] - r_i[2];
        // check if any atoms are at the singularity of the derivatives
        // occurs if x = y = 0 -> see discussion in Appendix B of J. Stimac's Dissertation
        // if so perturb all positions slightly
        if (i_coordinates[i][0] == 0.0 && i_coordinates[i][1] == 0.0) {
            rigid_shift = true;
        }
    }
    // Apply the rigid shift if necessary
    if (rigid_shift) {
        for (int i = 0; i < number_of_neighbors; i++) {
            i_coordinates[i][0] += 1e-6;
            i_coordinates[i][1] += 1e-6;
            i_coordinates[i][2] += 1e-6;
        }
    }

    //2. get spherical harmonics for center shifted coordinates
    auto Ylmi_real = std::vector<double>(number_of_neighbors * l_max_sq);
    auto Ylmi_imag = std::vector<double>(number_of_neighbors * l_max_sq);

    // 3. map coordinates to spherical coordinates
    auto i_coordinates_spherical = std::vector<double>(3 * MAX_NEIGHBORS);
    auto r_ij = std::vector<double>(MAX_NEIGHBORS);

    for (int i = 0; i < number_of_neighbors; i++) {
        r_ij[i] = std::sqrt(i_coordinates[i][0] * i_coordinates[i][0] + i_coordinates[i][1] * i_coordinates[i][1] +
                       i_coordinates[i][2] * i_coordinates[i][2]);
        Ylmi_all_l_from_r(q, r_ij.data(),
                              Ylmi_real.data() + i * l_max_sq,
                              Ylmi_imag.data() + i * l_max_sq);
    }

    // 4. get the radial basis function values
    auto gnl = std::vector<double>(number_of_neighbors * l_max_sq);
    bessel_basis(q, cutoff, number_of_neighbors, r_ij.data(), gnl.data());

    // 6. compute ct and st
    auto ct = std::vector<double>(l_max_sq * (2 * q + 1));
    auto st = std::vector<double>(l_max_sq * (2 * q + 1));
    for (int n = 0; n < q + 1; n++) {
        for (int l = 0; l < q + 1; l++) {
            for (int m = 0; m < 2 * q + 1; m++) {
                ct[n * (q + 1) * (2 * q + 1) + l * (2 * q + 1) + m] = 0.0;
                st[n * (q + 1) * (2 * q + 1) + l * (2 * q + 1) + m] = 0.0;
                for (int i = 0; i < number_of_neighbors; i++) {
                    ct[n * (q + 1) * (2 * q + 1) + l * (2 * q + 1) + m] +=
                            gnl[i * l_max_sq + n * (q + 1) + l] *
                            Ylmi_real[i * l_max_sq + m];
                    st[n * (q + 1) * (2 * q + 1) + l * (2 * q + 1) + m] +=
                            gnl[i * l_max_sq + n * (q + 1) + l] *
                            Ylmi_imag[i * l_max_sq + m];
                }
            }
        }
    }

    // 7. compute the descriptor
    for (int d = 0; d < this->width; d++) {
        int n0 = this->ln_params[d * 6 + 0];
        int n1 = this->ln_params[d * 6 + 1];
        int n2 = this->ln_params[d * 6 + 2];
        int l0 = this->ln_params[d * 6 + 3];
        int l1 = this->ln_params[d * 6 + 4];
        int l2 = this->ln_params[d * 6 + 5];
        double xi_desc = 0.0;
        for (int m0 = -l0; m0 < l0 + 1; m0++) {
            int m0_ind = m0 + q;
            int m1l = (m0 - l1 - l2 + std::abs(l1 - l2 + m0)) / 2;
            int m1n = (m0 + l1 + l2 - std::abs(l1 - l2 - m0)) / 2;
            for (int m1 = m1l; m1 < m1n + 1; m1++) {
                int m2 = m0 - m1;
                int m1_ind = m1 + q;
                int m2_ind = m2 + q;
                double C = clebsh_gordon(l1, l2, l0, m1, m2, m0);
                xi_desc += C * (ct[n0 * (q + 1) * (2 * q + 1) + l0 * (2 * q + 1) + m0_ind] *
                                (ct[n1 * (q + 1) * (2 * q + 1) + l1 * (2 * q + 1) + m1_ind] *
                                 ct[n2 * (q + 1) * (2 * q + 1) + l2 * (2 * q + 1) + m2_ind] -
                                 st[n1 * (q + 1) * (2 * q + 1) + l1 * (2 * q + 1) + m1_ind] *
                                 st[n2 * (q + 1) * (2 * q + 1) + l2 * (2 * q + 1) + m2_ind]) +
                                st[n0 * (q + 1) * (2 * q + 1) + l0 * (2 * q + 1) + m0_ind] *
                                (ct[n1 * (q + 1) * (2 * q + 1) + l1 * (2 * q + 1) + m1_ind] *
                                 st[n2 * (q + 1) * (2 * q + 1) + l2 * (2 * q + 1) + m2_ind] +
                                 st[n1 * (q + 1) * (2 * q + 1) + l1 * (2 * q + 1) + m1_ind] *
                                 ct[n2 * (q + 1) * (2 * q + 1) + l2 * (2 * q + 1) + m2_ind]));
            }
        }
        double coeff = std::pow(-1.0, l0) / std::sqrt(2 * l0 + 1);
        xi_desc *= coeff;
        desc[d] = xi_desc;
    }
}