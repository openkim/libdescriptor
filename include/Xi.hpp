#ifndef XI_HPP_
#define XI_HPP_

#include "Descriptors.hpp"
#include <iostream>

typedef double VectorOfSize3[3];

using namespace Descriptor;

class Xi final : public DescriptorKind{
public:
    Xi() {}; //TODO delete the default constructor
    Xi(std::string &filename);
    Xi(int q, double cutoff, std::vector<std::string> &species, std::string& radial_basis);

    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neighbor_lists,
                 int number_of_neighbors,
                 double *coordinates,
                 double *desc) override;

    void clone_empty(DescriptorKind *descriptorKind);

    void allocate_memory();

    int get_width();

    int q;
    double cutoff;
    std::vector<std::string> species_;
    std::string radial_basis = "bessel";
    ~Xi() override {} // all memory allocations are std::vector containers so no need to delete anything

private:
    std::vector<int> ln_params;
    std::vector<double> radial_basis_array;
    int l_max_sq;
    std::vector<double> clebsh_gordon_array;
    void create_clebsh_gordon_array();//TODO implement this function
    std::vector<double> Ylmi_real, Ylmi_imag;
    std::vector<double> center_shifted_neighbor_coordinates, center_shifted_neighbor_coordinates_zj;
    std::vector<double> i_coordinates_spherical, r_ij;
    std::vector<double> gnl;
    std::vector<double> ct, st;
};


#endif //XI_HPP_