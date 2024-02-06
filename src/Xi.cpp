#include "Xi.hpp"
#include "get_xi_parameters.hpp"
#include <cmath>
#include "maths/maths.hpp"
#include "file_io_utils.hpp"

#define MAX_NEIGHBORS 100

Xi::Xi(int q, double cutoff, std::vector<std::string> &species, std::string &radial_basis) {
    this->cutoff = cutoff;
    this->radial_basis = radial_basis;
    this->q = q;
    this->species_ = species;
    this->width = get_width();
    this->ln_params.resize(this->width * 6);
    copy_params(this->q, this->ln_params.data()); // get the parameters from the file }
    this->l_max_sq = (q + 1) * (q + 1);
    if (radial_basis != "bessel") {
        throw std::runtime_error("Radial basis " + radial_basis + " not implemented");
    }
    allocate_memory();
    //    this->create_clebsh_gordon_array();
}

void Xi::create_clebsh_gordon_array() {
}

int Xi::get_width() {
    return get_param_length(q);
}

void Xi::compute(int index,
                 int n_atoms,
                 int *species,
                 int *neighbor_lists,
                 int number_of_neighbors,
                 double *coordinates,
                 double *desc) {
    auto *coordinates_ = (VectorOfSize3 *) coordinates;
    std::array<double, 3> r_i = {coordinates_[index][0], coordinates_[index][1], coordinates_[index][2]};
    int ind_factor1 = (2 * q + 1);
    int ind_factor2 = (q + 1) * (2 * q + 1);

    //1. get center shifted coordinates for index atom. Neighbors exclude the center atom
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

    // 2. get spherical harmonics for center shifted coordinates
    // 3. map coordinates to spherical coordinates

    for (int i = 0; i < number_of_neighbors; i++) {
        r_ij[i] = std::sqrt(i_coordinates[i][0] * i_coordinates[i][0] + i_coordinates[i][1] * i_coordinates[i][1] +
                            i_coordinates[i][2] * i_coordinates[i][2]);
        std::array<double, 3> i_coords = {i_coordinates[i][0], i_coordinates[i][1], i_coordinates[i][2]};
        Ylmi_all_l_from_r(q, i_coords.data(),
                          Ylmi_real.data() + i * l_max_sq,
                          Ylmi_imag.data() + i * l_max_sq);
    }

    // 4. get the radial basis function values
    bessel_basis(q, cutoff, number_of_neighbors,
                 r_ij.data(), static_cast<int>(gnl.size()), gnl.data());

    // 5. compute ct and st
    for (int n = 0; n < q + 1; n++) {
        for (int l = 0; l < q + 1; l++) {
            for (int m = 0; m < 2 * l + 1; m++) {
                int idx = n * ind_factor2 + l * ind_factor1 + m;
                ct[idx] = 0.0;
                st[idx] = 0.0;
                for (int i = 0; i < number_of_neighbors; i++) {
                    ct[idx] +=
                            gnl[i * l_max_sq + n * (q + 1) + l] *
                            Ylmi_real[i * l_max_sq + l * l + m];
                    st[idx] +=
                            gnl[i * l_max_sq + n * (q + 1) + l] *
                            Ylmi_imag[i * l_max_sq + l * l + m];
                }
            }
        }
    }

    // 6. compute the descriptor
    for (int d = 0; d < width; d++) {
        int n0 = ln_params[d * 6 + 0];
        int n1 = ln_params[d * 6 + 1];
        int n2 = ln_params[d * 6 + 2];
        int l0 = ln_params[d * 6 + 3];
        int l1 = ln_params[d * 6 + 4];
        int l2 = ln_params[d * 6 + 5];
        double xi_desc = 0.0;
        for (int m0 = -l0; m0 < l0 + 1; m0++) {
            int m0_ind = m0 + l0;
            int m1l = (m0 - l1 - l2 + std::abs(l1 - l2 + m0)) / 2;
            int m1n = (m0 + l1 + l2 - std::abs(l1 - l2 - m0)) / 2;
            for (int m1 = m1l; m1 < m1n + 1; m1++) {
                int m2 = m0 - m1;
                int m1_ind = m1 + l1;
                int m2_ind = m2 + l2;
                double C = clebsh_gordon(l1, l2, l0, m1, m2, m0);
                xi_desc += C * (ct[n0 * ind_factor2 + l0 * ind_factor1 + m0_ind] *
                                (ct[n1 * ind_factor2 + l1 * ind_factor1 + m1_ind] *
                                 ct[n2 * ind_factor2 + l2 * ind_factor1 + m2_ind] -
                                 st[n1 * ind_factor2 + l1 * ind_factor1 + m1_ind] *
                                 st[n2 * ind_factor2 + l2 * ind_factor1 + m2_ind]) +
                                st[n0 * ind_factor2 + l0 * ind_factor1 + m0_ind] *
                                (ct[n1 * ind_factor2 + l1 * ind_factor1 + m1_ind] *
                                 st[n2 * ind_factor2 + l2 * ind_factor1 + m2_ind] +
                                 st[n1 * ind_factor2 + l1 * ind_factor1 + m1_ind] *
                                 ct[n2 * ind_factor2 + l2 * ind_factor1 + m2_ind]));
            }
        }
        double coeff = std::pow(-1.0, l0) / std::sqrt(2 * l0 + 1);
        xi_desc *= coeff;
        desc[d] = xi_desc;
    }
}

void Xi::clone_empty(DescriptorKind *descriptorKind) {
    auto *dxi = dynamic_cast<Xi *>(descriptorKind);
    q = dxi->q;
    cutoff = dxi->cutoff;
    species_ = dxi->species_;
    radial_basis = dxi->radial_basis;
    width = dxi->width;
    ln_params = dxi->ln_params;
    l_max_sq = dxi->l_max_sq;
    ln_params = dxi->ln_params;
    allocate_memory();
    // std::vector.resize sometimes does not work with Enzyme, so use assign instead
}

void Xi::allocate_memory() {
    Ylmi_real.resize(MAX_NEIGHBORS * (q + 1) * (q + 1));
    Ylmi_imag.resize(MAX_NEIGHBORS * (q + 1) * (q + 1));
    center_shifted_neighbor_coordinates.resize(3 * MAX_NEIGHBORS);
    center_shifted_neighbor_coordinates_zj.resize(3 * MAX_NEIGHBORS);
    i_coordinates_spherical.resize(3 * MAX_NEIGHBORS);
    r_ij.resize(MAX_NEIGHBORS);
    gnl.resize(MAX_NEIGHBORS * l_max_sq);
    ct.resize(l_max_sq * (2 * q + 1));
    st.resize(l_max_sq * (2 * q + 1));

}

Xi::Xi(std::string &filename) {

    // Open the descriptor file
    std::ifstream file = FileIOUtils::open_file(filename);

    // String containing data line and list of parameters
    std::vector<std::string> string_params;
    std::vector<double> double_params;
    std::vector<int> int_params;
    std::vector<bool> bool_params;
    std::string line;

    // params
    int q_, n_species_;
    double cutoff_;
    std::vector<std::string> species;
    std::string radial_basis_;

    // File format:
    /*
     * # q
     * 2
     * # cutoff
     * 5.0
     * # radial_basis
     * bessel
     * # species
     * 4
     * H C N O
     */

    // q
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_int_params(line, int_params, 1);
    q_ = int_params[0];
    int_params.clear();
    line.clear();

    // cutoff
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_double_params(line, double_params, 1);
    cutoff_ = double_params[0];
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
    this->cutoff = cutoff_;
    this->radial_basis = radial_basis_;
    this->q = q_;
    this->species_ = species;
    this->width = get_width();
    this->ln_params.resize(this->width * 6);
    copy_params(this->q, this->ln_params.data()); // get the parameters from the file }
    this->l_max_sq = (q + 1) * (q + 1);
    if (radial_basis != "bessel") {
        throw std::runtime_error("Radial basis " + radial_basis + " not implemented");
    }

    allocate_memory();
}
