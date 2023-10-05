#ifndef SOAP_HPP_
#define SOAP_HPP_

#include "helper.hpp"
#include "Descriptors.hpp"

#include <numeric>

#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include <memory>

#ifdef DIM
#undef DIM
#endif

#define DIM 3

typedef double VectorOfSizeDIM[DIM];

using namespace Descriptor;

class SOAP final : public DescriptorKind {
public:
    SOAP() {}; //TODO delete the default constructor
    SOAP(std::string &filename);

    SOAP(int n_max, int l_max, double cutoff, std::vector<std::string> &species, std::string radial_basis, double eta);

    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neighbor_lists,
                 int number_of_neighbors,
                 double *coordinates,
                 double *desc) override;

    void clone_empty(DescriptorKind *descriptorKind);

    void init_radial_basis_array();

    void allocate_memory();

    int get_width();

    int n_max;
    int l_max;
    double cutoff;
    double eta;
    int n_species;
    std::vector<std::string> species_;
    std::string radial_basis = "polynomial";
private:
    std::vector<double> radial_basis_array;
    int n_gl_quad_points = 100;
    std::vector<double> gl_quad_weights;
    std::vector<double> gl_quad_radial_grid_points;
    std::vector<double> gl_quad_radial_sq_grid_points;
    std::vector<double> Cznlm_real;
    std::vector<double> Cznlm_imag;
    std::vector<double> Cij_real;
    std::vector<double> Cij_imag;
    std::vector<double> power_spectrum;
    std::vector<double> Ylmi_real;
    std::vector<double> Ylmi_imag;
    std::vector<double> exp_eta_r2;
    std::vector<double> I_zj_real;
    std::vector<double> I_zj_imag;
    // std::vector<double> center_shifted_neighbor_coordinates_zj;
    // std::vector<double> center_shifted_neighbor_coordinates;

    int l_max_sq;
};

#endif // SOAP_HPP_