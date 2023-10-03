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

class SOAP final: public DescriptorKind {
public:
    SOAP() {}; //TODO delete the default constructor
    SOAP(std::string& filename);
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
    int width=-1;
    std::string radial_basis = "polynomial";
    std::vector<double> radial_basis_array;
    std::vector<double> gl_quad_weights;
    std::vector<double> gl_quad_radial_grid_points;
    std::vector<double> gl_quad_radial_sq_grid_points;
    int n_gl_quad_points = 100;
};

#endif // SOAP_HPP_