#include "SOAP.hpp"
#include "maths/maths.hpp"
#include "file_io_utils.hpp"
#include <memory>

#define MAX_NEIGHBORS 100

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
    this->l_max_sq = (l_max + 1) * (l_max + 1);
    allocate_memory();
    init_radial_basis_array();
    this->width = get_width();
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
        polynomial_basis(n_max, cutoff, n_gl_quad_points, gl_quad_radial_grid_points.data(),
                         radial_basis_array.data());
        for (int i = 0; i < n_gl_quad_points; i++) {
            gl_quad_radial_sq_grid_points[i] = gl_quad_radial_grid_points[i] * gl_quad_radial_grid_points[i];
            exp_eta_r2[i] = std::exp(-1 * eta * gl_quad_radial_sq_grid_points[i]);
        }
    } else {
        throw std::invalid_argument("radial_basis must be one of: polynomial");
    }
}

void
SOAP::compute(int index, int n_atoms, int *species, int *neighbor_lists, int number_of_neighbors, double *coordinates,
              double *desc) {
    auto *coordinates_ = (VectorOfSizeDIM *) coordinates;
    std::array<double, 3> r_i = {coordinates_[index][0], coordinates_[index][1], coordinates_[index][2]};

    // TODO: these 2 coord allocations are needed in every call, else Enzyme fails/gives different results
    // Restest them with later versions of Enzyme
    std::fill(I_zj_real.begin(), I_zj_real.end(), 0.0);
    std::fill(I_zj_imag.begin(), I_zj_imag.end(), 0.0);
    auto center_shifted_neighbor_coordinates = std::vector<double>(3 * MAX_NEIGHBORS);
    auto center_shifted_neighbor_coordinates_zj = std::vector<double>(3 * MAX_NEIGHBORS);


    //1. get center shifted coordinates for index atom. Neighors exclude the center atom
    auto i_coordinates = (VectorOfSizeDIM *) center_shifted_neighbor_coordinates.data();
    for (int i = 0; i < number_of_neighbors; i++) {
        int neighbor_index = neighbor_lists[i];
        std::array<double, 3> r_j = {coordinates_[neighbor_index][0], coordinates_[neighbor_index][1],
                                     coordinates_[neighbor_index][2]};
        i_coordinates[i][0] = r_j[0] - r_i[0];
        i_coordinates[i][1] = r_j[1] - r_i[1];
        i_coordinates[i][2] = r_j[2] - r_i[2];
    }
    //2. get all the neighbor coordinates of a given species
    for (int zj = 0; zj < n_species; zj++) {
        int n_neighbors_zj = 0;
        std::vector<int> species_j_indices(MAX_NEIGHBORS, -1);
        int i_neigh_j = 0;
        for (int i = 0; i < number_of_neighbors; i++) {
            if (species[neighbor_lists[i]] == zj) {
                species_j_indices[i_neigh_j] = i;
                n_neighbors_zj++;
                i_neigh_j++;
            }
        }
        auto j_coordinates = (VectorOfSizeDIM *) center_shifted_neighbor_coordinates_zj.data();
        for (int j = 0; j < n_neighbors_zj; j++) {
            j_coordinates[j][0] = i_coordinates[species_j_indices[j]][0];
            j_coordinates[j][1] = i_coordinates[species_j_indices[j]][1];
            j_coordinates[j][2] = i_coordinates[species_j_indices[j]][2];
        }

        //3. For each species element, get Ylmi and store
        //4. For each species element, get I_ij and store
        std::fill(Ylmi_real.begin(), Ylmi_real.begin() + l_max_sq * n_neighbors_zj, 0.0);
        std::fill(Ylmi_imag.begin(), Ylmi_imag.begin() + l_max_sq * n_neighbors_zj, 0.0);
        for (int j = 0; j < n_neighbors_zj; j++) {
            std::array<double, 3> r_j = {j_coordinates[j][0], j_coordinates[j][1], j_coordinates[j][2]};

            Ylmi_all_l_from_r(l_max, r_j.data(),
                              Ylmi_real.data() + j * l_max_sq,
                              Ylmi_imag.data() + j * l_max_sq);

            double r_ij_sq = r_j[0] * r_j[0] + r_j[1] * r_j[1] + r_j[2] * r_j[2];
            double r_ij = std::sqrt(r_ij_sq);

            //5. Sum over all the I_ij of a given species element
            double two_eta_r2_rij = 0;
            int lm_idx = -1;
            for (int l = 0; l < l_max + 1; l++) {
                for (int m = -l; m <= l; m++) {
                    lm_idx = l * l + l + m;
                    for (int i = 0; i < n_gl_quad_points; i++) {
                        int idx = n_gl_quad_points * lm_idx + i;
                        two_eta_r2_rij = spherical_in(l, 2 * eta * gl_quad_radial_grid_points[i] * r_ij) *
                                         std::exp(-1 * eta * r_ij_sq) * exp_eta_r2[i];
                        I_zj_real[idx] += two_eta_r2_rij * Ylmi_real[j * l_max_sq + lm_idx];
                        I_zj_imag[idx] += two_eta_r2_rij * Ylmi_imag[j * l_max_sq + lm_idx];
                    }
                }
            }
        }

        //6. numerically integrate over the radial grid points
        //Cij_real = \int r^2 R(r) * I_zj_real, Cij_imag = \int r^2 R(r) * I_zj_imag = \sum_i r[i]^2 R(r[i]) * I_zj_imag[i] w_i
        std::fill(Cij_real.begin(), Cij_real.end(), 0.0);
        std::fill(Cij_imag.begin(), Cij_imag.end(), 0.0);
        int lm_idx = -1;
        double gl_weigth_r2 = 0;
        for (int n = 0; n < n_max; n++) {
            for (int l = 0; l < l_max + 1; l++) {
                for (int m = -l; m <= l; m++) {
                    lm_idx = l * l + l + m;
                    for (int i = 0; i < n_gl_quad_points; i++) {
                        gl_weigth_r2 = radial_basis_array[n * n_gl_quad_points + i] * gl_quad_radial_sq_grid_points[i] *
                                       gl_quad_weights[i];
                        Cij_real[n * l_max_sq + lm_idx] +=
                                gl_weigth_r2 * I_zj_real[l * (l + 1) * n_gl_quad_points + m * n_gl_quad_points + i];
                        Cij_imag[n * l_max_sq + lm_idx] +=
                                gl_weigth_r2 * I_zj_imag[l * (l + 1) * n_gl_quad_points + m * n_gl_quad_points + i];
                    }
                    //7. Store the Cznlm for each species element
                    Cznlm_real[zj * n_max * l_max_sq + n * l_max_sq + lm_idx] = Cij_real[n * l_max_sq + lm_idx];
                    Cznlm_imag[zj * n_max * l_max_sq + n * l_max_sq + lm_idx] = Cij_imag[n * l_max_sq + lm_idx];
                }
            }
        }
    }

    //8. compute the power spectrum for each species element
    int tmp_index = 0;//TODO: get proper index
    int lm_idx = -1;
    std::fill(power_spectrum.begin(), power_spectrum.end(), 0.0);
    for (int z_j = 0; z_j < n_species; z_j++) {
        for (int z_i = z_j; z_i < n_species; z_i++) {
            for (int n_j = 0; n_j < n_max; n_j++) {
                for (int n_i = n_j; n_i < n_max; n_i++) {
                    for (int l = 0; l < l_max + 1; l++) {
                        for (int m = -l; m <= l; m++) {
                            lm_idx = l * l + l + m;
                            power_spectrum[tmp_index] +=
                                    Cznlm_real[z_j * n_max * l_max_sq + n_j * l_max_sq + lm_idx] *
                                    Cznlm_real[z_i * n_max * l_max_sq + n_i * l_max_sq + lm_idx] +
                                    Cznlm_imag[z_j * n_max * l_max_sq + n_j * l_max_sq + lm_idx] *
                                    Cznlm_imag[z_i * n_max * l_max_sq + n_i * l_max_sq + lm_idx];
                        }
                        power_spectrum[tmp_index] *= 124.02510672119926 * std::sqrt(
                                8.0 / (2.0 * l + 1.0)); //  M_PI * 4 * std::pow(M_PI, 2) = 124.02510672119926
                        desc[tmp_index] = power_spectrum[tmp_index];
                        tmp_index++;
                    }
                }
            }
        }
    }
}


void SOAP::allocate_memory() {

    //1. spherical bessel I_ij functions (l_max + 1)**2 * n_gl_quad_points  (real and imaginary part)
    I_zj_real = std::vector<double>(l_max_sq * n_gl_quad_points); // preallocated
    I_zj_imag = std::vector<double>(l_max_sq * n_gl_quad_points); // preallocated

    //2. Projection Cnlm coefficients
    Cznlm_real = std::vector<double>(n_species * n_max * l_max_sq);
    Cznlm_imag = std::vector<double>(n_species * n_max * l_max_sq);
    Cij_real = std::vector<double>(n_max * l_max_sq);
    Cij_imag = std::vector<double>(n_max * l_max_sq);

    //3. spherical harmonics: first (l_max + 1)**2 is real part, second (l_max + 1)**2 is imaginary part
    Ylmi_real = std::vector<double>(l_max_sq * MAX_NEIGHBORS);
    Ylmi_imag = std::vector<double>(l_max_sq * MAX_NEIGHBORS);

    //4. exp(-eta * r**2) for each grid point for integration
    exp_eta_r2 = std::vector<double>(n_gl_quad_points);
    gl_quad_radial_sq_grid_points = std::vector<double>(n_gl_quad_points);

    // center_shifted_neighbor_coordinates = std::vector<double>(3 * MAX_NEIGHBORS);
    // center_shifted_neighbor_coordinates_zj = std::vector<double>(3 * MAX_NEIGHBORS);
    //5. power spectrum
    power_spectrum = std::vector<double>(n_species * (n_species + 1) / 2 * n_max * (n_max + 1) / 2 * (l_max + 1));

}

void SOAP::clone_empty(DescriptorKind *descriptorKind) {
    auto d_soap = dynamic_cast<SOAP *>(descriptorKind);
    n_max = d_soap->n_max;
    l_max = d_soap->l_max;
    cutoff = d_soap->cutoff;
    n_species = d_soap->n_species;
    eta = d_soap->eta;
    l_max_sq = d_soap->l_max_sq;
    allocate_memory();
    init_radial_basis_array();
    width = d_soap->width;
    // zero out the memory
    //    if (radial_basis == "polynomial") {
    //        for (int i = 0; i < d_soap->radial_basis_array.size(); i++) {
    //            radial_basis_array[i] = 0.0;
    //            gl_quad_radial_sq_grid_points[i] = 0.0;
    //            gl_quad_weights[i] = 0.0;
    //            gl_quad_radial_grid_points[i] = 0.0;
    //        }
    //    } else {
    //        throw std::invalid_argument("radial_basis must be one of: polynomial");
    //    }
}

SOAP::SOAP(std::string &filename) {

    // Open the descriptor file
    std::ifstream file = FileIOUtils::open_file(filename);

    // String containing data line and list of parameters
    std::vector<std::string> string_params;
    std::vector<double> double_params;
    std::vector<int> int_params;
    std::vector<bool> bool_params;
    std::string line;

    // params
    int n_max_, l_max_, n_species_;
    double cutoff_, eta_;
    std::vector<std::string> species;
    std::string radial_basis_;

    // File format:
    /*
     * # n_max
     * 2
     * # l_max
     * 2
     * # cutoff
     * 5.0
     * # eta
     * 1.0
     * # radial_basis
     * polynomial
     * # species
     * 4
     * H C N O
     */
    // n_max
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_int_params(line, int_params, 1);
    n_max_ = int_params[0];
    int_params.clear();
    line.clear();

    // l_max
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_int_params(line, int_params, 1);
    l_max_ = int_params[0];
    int_params.clear();
    line.clear();

    // cutoff
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_double_params(line, double_params, 1);
    cutoff_ = double_params[0];
    double_params.clear();
    line.clear();

    // eta
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_double_params(line, double_params, 1);
    eta_ = double_params[0];
    double_params.clear();
    line.clear();

    // radial_basis
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_string_params(line, string_params, 1);
    radial_basis_ = string_params[0];
    string_params.clear();
    line.clear();

    // species
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_int_params(line, int_params, 1);
    n_species_ = int_params[0];
    int_params.clear();
    line.clear();

    FileIOUtils::parse_string_params(line, string_params, n_species_);
    for (auto &species_name: string_params) {
        species.push_back(species_name);
    }
    string_params.clear();
    line.clear();

    file.close();

    // Initialize the descriptor
    this->n_max = n_max_;
    this->l_max = l_max_;
    this->cutoff = cutoff_;
    this->species_ = species;
    this->n_species = species.size();
    this->radial_basis = std::move(radial_basis);
    this->eta = eta_;
    l_max_sq = (l_max + 1) * (l_max + 1);
    init_radial_basis_array();
    allocate_memory();
    this->width = this->get_width();

}

