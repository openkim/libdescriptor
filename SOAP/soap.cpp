#include "soap.hpp"
#include "maths/maths.hpp"
#include <memory>

using namespace Descriptor;

SOAP::SOAP(int n_max, int l_max, double cutoff, std::vector<std::string> &species, std::string radial_basis,
           double eta) {
    this->n_max = n_max;
    this->l_max = l_max;
    this->cutoff = cutoff;
    this->species_ = species;
    this->n_species = species.size();
    this->radial_basis = std::move(radial_basis);
    this->eta = eta;
    init_radial_basis_array();
    allocate_memory();

}

int SOAP::get_width() {
    if (width == -1) {
        width = (((n_species + 1) * n_species) / 2 * (n_max * (n_max + 1)) * (l_max + 1)) / 2;
    }
    return this->width;
}

void SOAP::init_radial_basis_array() {
    if (radial_basis == "polynomial") {
        gl_quad_weights = get_gl_weights();
        gl_quad_radial_grid_points = get_gl_grid(cutoff);
        n_gl_quad_points = gl_quad_weights.size();

        radial_basis_array = std::vector<double>(n_max * n_gl_quad_points);
        polynomial_basis(n_max, cutoff, n_gl_quad_points, gl_quad_radial_grid_points.data(), radial_basis_array.data());
        gl_quad_radial_sq_grid_points = std::vector<double>(n_gl_quad_points);
    } else {
        throw std::invalid_argument("radial_basis must be one of: polynomial");
    }
}

// Fix this later
//void SOAP::compute(int index,
//                   int n_atoms,
//                   int *species,
//                   int *neighbor_lists,
//                   int number_of_neighbors,
//                   double *coordinates,
//                   double *desc) {
//    auto *coordinates_ = (VectorOfSizeDIM *) coordinates;
//    std::array<double, 3> r_i = {coordinates_[index][0], coordinates_[index][1], coordinates_[index][2]};
//    int species_i = species[index];
//
//    // All preallocations of arrays
//    //1. spherical harmonics: first (l_max + 1)**2 is real part, second (l_max + 1)**2 is imaginary part
//    std::vector<double> Ylmi(2 * (l_max + 1) * (l_max + 1)); // preallocate and init to zero
//    // std::vector<double> Ylmi_real((l_max + 1) * (l_max + 1)), Ylmi_imag((l_max + 1) * (l_max + 1));
//    //2. exp(-eta * r**2) for each grid point for integration
//    std::vector<double> exp_radial_grid_factors(n_gl_quad_points);
//    //3. spherical bessel factors for each grid point for integration i(2arr_i)
//    std::vector<double> spherical_bessel_factors((l_max + 1) * n_gl_quad_points); // preallocated
//    //4. I_ij functions (l_max + 1)**2 * n_gl_quad_points  (real and imaginary part)
//    std::vector<double> I_ij_real((l_max + 1) * (l_max + 1) * n_gl_quad_points); // preallocated
//    std::vector<double> I_ij_imag((l_max + 1) * (l_max + 1) * n_gl_quad_points); // preallocated
//    //5. Cnlm and Cnlm_conj for each neighbor
//    std::vector<double> Cznlm_real((l_max + 1) * (l_max + 1) * n_max * n_species); // preallocated
//    std::vector<double> Cznlm_imag((l_max + 1) * (l_max + 1) * n_max * n_species); // preallocated
//    //6. power_spectrum of size
//    std::vector<double> power_spectrum((n_max * n_species) * (n_max * n_species + 1) / 2 * (l_max + 1))
//); // preallocated
//
//    std::fill(Ylmi.begin(), Ylmi.end(), 0.0);
//    int Ylmi_size = (l_max + 1) * (l_max + 1);
//
//    for(int i = 0; i < n_gl_quad_points; i++){
//        double r = gl_quad_radial_grid_points[i];
//        double r2 = r * r;
//        gl_quad_radial_sq_grid_points[i] = r2;
//        exp_radial_grid_factors[i] = std::exp(-eta * r2);
//    }
//
//    for (int i = 0; i < number_of_neighbors; i++) {
//        int j = neighbor_lists[i];
//        int species_j = species[j];
//        std::array<double, 3> r_j = {coordinates_[j][0], coordinates_[j][1], coordinates_[j][2]};
//        std::array<double, 3> r_ij = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};
//        double r_ij_norm_sq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
//        double r_ij_norm = std::sqrt(r_ij_norm_sq);
//        Ylmi_all_l_from_r(l_max, r_ij.data(), Ylmi.data());
//        double exp_rij2 = std::exp(-eta * r_ij_norm_sq);
//        for (int l = 0; l < l_max + 1; l++){
//            for(int k = 0; k < n_gl_quad_points; k++){
//                spherical_bessel_factors[l * n_gl_quad_points + k] =
//                        spherical_in(l, 2 * eta * r_ij_norm * gl_quad_radial_grid_points[k]);
//            }
//        }
//        for (int l = 0; l < l_max + 1; l++){
//            int l_index = l * (l + 1);
//            for (int m = -l; m < l + 1; m++){
//                int lm = l_index + m;
//                double Ylmi_real = Ylmi[lm];
//                double Ylmi_imag = Ylmi[lm + Ylmi_size];
//                for(int k = 0; k < n_gl_quad_points; k++){
//                    double spherical_bessel_factor = spherical_bessel_factors[l * n_gl_quad_points + k];
//                    double I_ij_real_ = spherical_bessel_factor * Ylmi_real;
//                    double I_ij_imag_ = spherical_bessel_factor * Ylmi_imag;
//                    I_ij_real[l * n_gl_quad_points + k] = I_ij_real_;
//                    I_ij_imag[l * n_gl_quad_points + k] = I_ij_imag_;
//                }
//            }
//        }
//        // integration over radial grid points
//        for (int l = 0; l < l_max + 1; l++) {
//            int l_index = l * (l + 1);
//            for (int m = -l; m < l + 1; m++) {
//                int lm = l_index + m;
//                double Cnlm_real_ = 0.0;
//                double Cnlm_imag_ = 0.0;
//                for (int n = 0; n < n_max; n++) {
//                    for (int k = 0; k < n_gl_quad_points; k++) {
//                        double I_ij_real_ = I_ij_real[l * n_gl_quad_points + k];
//                        double I_ij_imag_ = I_ij_imag[l * n_gl_quad_points + k];
//                        double exp_r2 = exp_radial_grid_factors[k];
//                        // C = sum_k i(2 a r rij) * Ylmi(rij) * exp(-eta * r**2) * exp(-eta * r_ij**2) * g_n(r) * r**2
//                        Cnlm_real_ += gl_quad_weights[k] * I_ij_real_ * exp_r2 * exp_rij2 *
//                                      radial_basis_array[n * n_gl_quad_points + k] * gl_quad_radial_sq_grid_points[k];
//                        Cnlm_imag_ += gl_quad_weights[k] * I_ij_imag_ * exp_r2 * exp_rij2 *
//                                      radial_basis_array[k * n_gl_quad_points + k] * gl_quad_radial_sq_grid_points[k];
//                    }
//                    index = species_j * n_max * (l_max + 1) * (l_max + 1) +
//                            n * (l_max + 1) * (l_max + 1) +
//                            l * (l + 1) + m;
//                    Cznlm_real[index] = Cnlm_real_;
//                    Cznlm_imag[index] = Cnlm_imag_;
//                }
//            }
//        }
//    }
//    int p_index = 0;
//    for (int z1 = 0; z1 < n_species; z1++){
//        for (int z2 = 0; z2 < n_species; z2++){
//            for (int n1 = 0; n1 < n_max; n1++){
//                for (int n2 = 0; n2 < n_max; n2++){
//                    for (int l = 0; l < l_max + 1; l++){
//                        int l_index = l * (l + 1);
//                        for (int m = -l; m < l + 1; m++){
//                            if ((z1 >= z2) && (n1 >= n2)){
//                                power_spectrum[p_index] = Cznlm_real[z1 * n_max * (l_max + 1) * (l_max + 1) +
//                                                                      n1 * (l_max + 1) * (l_max + 1) +
//                                                                      l * (l + 1) + m] *
//                                                          Cznlm_real[z2 * n_max * (l_max + 1) * (l_max + 1) +
//                                                                      n2 * (l_max + 1) * (l_max + 1) +
//                                                                      l * (l + 1) + m] +
//                                                          Cznlm_imag[z1 * n_max * (l_max + 1) * (l_max + 1) +
//                                                                      n1 * (l_max + 1) * (l_max + 1) +
//                                                                      l * (l + 1) + m] *
//                                                          Cznlm_imag[z2 * n_max * (l_max + 1) * (l_max + 1) +
//                                                                      n2 * (l_max + 1) * (l_max + 1) +
//                                                                      l * (l + 1) + m];
//                                p_index++;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//}

void
SOAP::compute(int index, int n_atoms, int *species, int *neighbor_lists, int number_of_neighbors, double *coordinates,
              double *desc) {
    auto *coordinates_ = (VectorOfSizeDIM *) coordinates;
    std::array<double, 3> r_i = {coordinates_[index][0], coordinates_[index][1], coordinates_[index][2]};

    std::vector<double> Cznlm_real = std::vector<double>(n_species * n_max * (l_max + 1) * (l_max + 1));
    std::vector<double> Cznlm_imag = std::vector<double>(n_species * n_max * (l_max + 1) * (l_max + 1));

    //1. get center shifted coordinates for index atom. Neighors exclude the center atom
    double *center_shifted_neighbor_coordinates = new double[3 * number_of_neighbors];
    auto i_coordinates = (VectorOfSizeDIM *) center_shifted_neighbor_coordinates;
    for (int i = 0; i < number_of_neighbors; i++) {
        int neighbor_index = neighbor_lists[i];
        std::array<double, 3> r_j = {coordinates_[neighbor_index][0], coordinates_[neighbor_index][1],
                                     coordinates_[neighbor_index][2]};
        i_coordinates[i][0] = r_j[0] - r_i[0];
        i_coordinates[i][1] = r_j[0] - r_i[1];
        i_coordinates[i][2] = r_j[0] - r_i[2];
    }
    //2. get all the neighbor coordinates of a given species
    int zi = species[index];
    for (int zj = 0; zj < n_species; zj++) {
        int n_neighbors_zj = 0;
        std::vector<int> species_j_indices;
        for (int i = 0; i < number_of_neighbors; i++) {
            if (species[neighbor_lists[i]] == zj) {
                species_j_indices.push_back(i);
                n_neighbors_zj++;
            }
        }
        double *center_shifted_neighbor_coordinates_zj = new double[3 * n_neighbors_zj];
        auto j_coordinates = (VectorOfSizeDIM *) center_shifted_neighbor_coordinates_zj;
        for (int j = 0; j < n_neighbors_zj; j++) {
            j_coordinates[j][0] = i_coordinates[species_j_indices[j]][0];
            j_coordinates[j][1] = i_coordinates[species_j_indices[j]][1];
            j_coordinates[j][2] = i_coordinates[species_j_indices[j]][2];
        }

        //3. For each species element, get Ylmi and store
        //4. For each species element, get I_ij and store
        std::vector<double> Ylmi_real(
                (l_max + 1) * (l_max + 1) * n_neighbors_zj); // preallocate and init to zero, 2 for imag and real
        std::vector<double> Ylmi_imag(
                (l_max + 1) * (l_max + 1) * n_neighbors_zj); // preallocate and init to zero, 2 for imag and real
        std::vector<double> bessel_i_zj(n_neighbors_zj * (l_max + 1) *
                                        n_gl_quad_points); // modeified bessel function, i(l, 2arri) for each ri each r
        std::fill(Ylmi_real.begin(), Ylmi_real.end(), 0.0);
        std::fill(Ylmi_imag.begin(), Ylmi_imag.end(), 0.0);
        for (int j = 0; j < n_neighbors_zj; j++) {
            std::array<double, 3> r_j = {j_coordinates[j][0], j_coordinates[j][1], j_coordinates[j][2]};

            Ylmi_all_l_from_r(l_max, r_j.data(),
                              Ylmi_real.data() + j * (l_max + 1) * (l_max + 1),
                              Ylmi_imag.data() + j * (l_max + 1) * (l_max +
                                                                    1)); // 2 for imag and real, org: real_1,imag_1,real_2,imag_2,real_3,imag_3 ..

            double r_ij = std::sqrt(r_j[0] * r_j[0] + r_j[1] * r_j[1] + r_j[2] * r_j[2]);
            // bessel_i_zy
            for (int l = 0; l < l_max + 1; l++) {
                for (int i = 0; i < n_gl_quad_points; i++) {
                    double r_sq = gl_quad_radial_grid_points[i] * gl_quad_radial_grid_points[i];
                    bessel_i_zj[j * (l_max + 1) * n_gl_quad_points + l * n_gl_quad_points + i] = spherical_in(l,
                                                                                                              2 * eta *
                                                                                                              r_sq *
                                                                                                              r_ij);
                }
            }
        }

        //5. Sum over all the I_ij of a given species element
        std::vector<double> I_zj_real(
                (l_max + 1) * (l_max + 1) * n_gl_quad_points); // Summed i(l, 2arri) exp(-arij^2) exp(-ar^2) Ylm(ri)
        std::vector<double> I_zj_imag((l_max + 1) * (l_max + 1) * n_gl_quad_points);
        for (int j = 0; j < n_neighbors_zj; j++) {
            double r_ij = std::sqrt(
                    j_coordinates[j][0] * j_coordinates[j][0] +
                    j_coordinates[j][1] * j_coordinates[j][1] +
                    j_coordinates[j][2] * j_coordinates[j][2]);
            for (int l = 0; l < l_max + 1; l++) {
                for (int m = -l; m <= l; m++) {
                    for (int i = 0; i < n_gl_quad_points; i++) {
                        I_zj_real[n_gl_quad_points * (l * (l + 1) + m) + i] +=
                                bessel_i_zj[j * (l_max + 1) * n_gl_quad_points + l * n_gl_quad_points + i] *
                                std::exp(-1 * eta * r_ij * r_ij) *
                                Ylmi_real[j * (l_max + 1) * (l_max + 1) + l * l + l + m] *
                                std::exp(-1 * eta * gl_quad_radial_grid_points[i] *
                                         gl_quad_radial_grid_points[i]);
                        I_zj_imag[n_gl_quad_points * (l * (l + 1) + m) + i] +=
                                bessel_i_zj[j * (l_max + 1) * n_gl_quad_points + l * n_gl_quad_points + i] *
                                std::exp(-1 * eta * r_ij * r_ij) *
                                Ylmi_imag[j * (l_max + 1) * (l_max + 1) + l * l + l + m] *
                                std::exp(-1 * eta * gl_quad_radial_grid_points[i] *
                                         gl_quad_radial_grid_points[i]);
                    }
                }
            }
        }

        //6. numerically integrate over the radial grid points
        // Cij_real = \int r^2 R(r) * I_zj_real, Cij_conj = \int r^2 R(r) * I_zj_imag = \sum_i r[i]^2 R(r[i]) * I_zj_imag[i] w_i
        std::vector<double> Cij_real(n_max * (l_max + 1) * (l_max + 1));
        std::vector<double> Cij_conj(n_max * (l_max + 1) * (l_max + 1));
        std::fill(Cij_real.begin(), Cij_real.end(), 0.0);
        std::fill(Cij_conj.begin(), Cij_conj.end(), 0.0);
        for (int n = 0; n < n_max; n++) {
            for (int l = 0; l < l_max + 1; l++) {
                for (int m = -l; m <= l; m++) {
                    for (int i = 0; i < n_gl_quad_points; i++) {
                        Cij_real[n * (l_max + 1) * (l_max + 1) + l * (l + 1) + m] +=
                                gl_quad_radial_grid_points[i] * gl_quad_radial_grid_points[i] *
                                radial_basis_array[n * n_gl_quad_points + i] *
                                I_zj_real[l * (l + 1) * n_gl_quad_points + m * n_gl_quad_points + i] *
                                gl_quad_weights[i];
                        Cij_conj[n * (l_max + 1) * (l_max + 1) + l * (l + 1) + m] +=
                                gl_quad_radial_grid_points[i] * gl_quad_radial_grid_points[i] *
                                radial_basis_array[n * n_gl_quad_points + i] *
                                I_zj_imag[l * (l + 1) * n_gl_quad_points + m * n_gl_quad_points + i] *
                                gl_quad_weights[i];
                    }
                    //7. Store the Cznlm for each species element
                    Cznlm_real[zj * n_max * (l_max + 1) * (l_max + 1) + n * (l_max + 1) * (l_max + 1) + l * (l + 1) +
                               m] = Cij_real[n * (l_max + 1) * (l_max + 1) + l * (l + 1) + m];
                    Cznlm_imag[zj * n_max * (l_max + 1) * (l_max + 1) + n * (l_max + 1) * (l_max + 1) + l * (l + 1) +
                               m] = Cij_conj[n * (l_max + 1) * (l_max + 1) + l * (l + 1) + m];
                }
            }
        }
        delete[] center_shifted_neighbor_coordinates_zj;
    }

    //8. compute the power spectrum for each species element
    std::vector<double> power_spectrum(n_species * (n_species + 1) / 2 * n_max * (n_max + 1) / 2 * (l_max + 1));
    int tmp_index = 0;//TODO: get proper index
    std::fill(power_spectrum.begin(), power_spectrum.end(), 0.0);
    for (int z_j = 0; z_j < n_species; z_j++) {
        for (int z_i = z_j; z_i < n_species; z_i++) {
            for (int n_j = 0; n_j < n_max; n_j++) {
                for (int n_i = n_j; n_i < n_max; n_i++) {
                    for (int l = 0; l < l_max + 1; l++) {
                        for (int m = -l; m <= l; m++) {
                            power_spectrum[tmp_index] +=
                                    Cznlm_real[z_j * n_max * (l_max + 1) * (l_max + 1) +
                                               n_j * (l_max + 1) * (l_max + 1) + l * (l + 1) + m] *
                                    Cznlm_real[z_i * n_max * (l_max + 1) * (l_max + 1) +
                                               n_i * (l_max + 1) * (l_max + 1) + l * (l + 1) + m] +
                                    Cznlm_imag[z_j * n_max * (l_max + 1) * (l_max + 1) +
                                               n_j * (l_max + 1) * (l_max + 1) + l * (l + 1) + m] *
                                    Cznlm_imag[z_i * n_max * (l_max + 1) * (l_max + 1) +
                                               n_i * (l_max + 1) * (l_max + 1) + l * (l + 1) + m];
                        }
                        power_spectrum[tmp_index] *= M_PI * std::sqrt(8.0 / (2.0 * l + 1.0)) * 4 * std::pow(M_PI, 2);
                        desc[tmp_index] = power_spectrum[tmp_index];
                        tmp_index++;
                    }
                }
            }
        }
    }

    delete[] center_shifted_neighbor_coordinates;
}


void SOAP::allocate_memory() {
    //1. spherical harmonics: first (l_max + 1)**2 is real part, second (l_max + 1)**2 is imaginary part
    //    Ylmi = std::vector<double>(2 * (l_max + 1) * (l_max + 1)); // preallocate and init to zero
    //    Ylmi_real = std::vector<double>((l_max + 1) * (l_max + 1));
    //    Ylmi_imag = std::vector<double>((l_max + 1) * (l_max + 1));
    //    //2. exp(-eta * r**2) for each grid point for integration
    //    exp_radial_grid_factors = std::vector<double>(n_gl_quad_points);
    //    //3. spherical bessel factors for each grid point for integration i(2arr_i)
    //    spherical_bessel_factors = std::vector<double>((l_max + 1) * n_gl_quad_points); // preallocated
    //    //4. I_ij functions (l_max + 1)**2 * n_gl_quad_points  (real and imaginary part)
    //    I_ij_real = std::vector<double>((l_max + 1) * (l_max + 1) * n_gl_quad_points); // preallocated
    //    I_ij_imag = std::vector<double>((l_max + 1) * (l_max + 1) * n_gl_quad_points); // preallocated
}