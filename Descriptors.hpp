// Generic header for inclusion and library design

#ifndef __DESCRIPTOR_HPP_
#define __DESCRIPTOR_HPP_

#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

typedef autodiff::VectorXreal double_vector;
typedef autodiff::MatrixXreal double_matrix;
typedef autodiff::real double_scalar;

namespace Descriptor {
    enum AvailableDescriptor {
        KindSymmetryFunctions, //!< For selecting Behler Symmetry Functions (SymmetryFunctions)
        KindBispectrum //!< For selecting Bispectrum descriptor
    };

    class DescriptorKind;
    inline void fwd_diff(int n_atoms, int *species, int *neighbor_list, int *number_of_neighs,
                         double *coordinates, double *desc,
                         DescriptorKind *descriptor_kind) { throw std::logic_error("Function not implemented yet."); }
    void gradient(int n_contributing_atoms, int n_total_atoms,  int *species, int *neighbor_list, int *number_of_neighs,
                  double *coordinates, double *d_coordinates, double *desc,
                  double *dE_dzeta, DescriptorKind *descriptor_to_diff);
    void gradient_single_atom(int index, int n_total_atoms, int *species, int *neighbor_list, int number_of_neighs,
                              double *coordinates, double *d_coordinates, double *desc,
                              double *dE_dzeta, DescriptorKind *descriptor_to_diff);

    void num_gradient_single_atom(int index, int n_contributing_atoms, int n_total_atoms , int *species, int *neighbor_list, int number_of_neighs,
                              double *coordinates, double *d_coordinates, double *dE_dzeta,
                              DescriptorKind *descriptor_to_diff);
    inline void jacobian(int n_atoms, int *species, int *neighbor_list, int *number_of_neighs,
                         double *coordinates, double *d_coordinates, double *desc,
                         double *dzeta_dr, DescriptorKind *descriptor_to_diff) {

        throw std::logic_error("Function not implemented yet");
    }

    void compute(int n_contributing_atoms, int n_total_atoms , int *species, int *neighbor_list, int *number_of_neighs,
                 double *coordinates, double *desc,
                 DescriptorKind *descriptor_kind);

    void compute_single_atom(int index, int n_total_atoms, int *species, int *number_of_neighs, int number_of_neigh_list,
                             double *coordinates, double *desc,
                             DescriptorKind *descriptor_kind);
}

class Descriptor::DescriptorKind {

public:
    AvailableDescriptor descriptor_kind; //!< Kind of instantiated descriptor, will be used in creating clone for AD
    std::string descriptor_param_file; //!< Full path to descriptor parameter file.
    int width=-1; //!< Dimension of the descriptor

    DescriptorKind() = default;

    static DescriptorKind *
    initDescriptor(AvailableDescriptor availableDescriptorKind); //!< Initialize an empty descriptor of a kind.
    static DescriptorKind *initDescriptor(std::string &file_name,
                                          AvailableDescriptor availableDescriptorKind);
    virtual void compute(int index,
                         int n_atoms,
                         int *species,
                         int *neighbor_lists,
                         int number_of_neighbors,
                         double_vector &coordinates,
                         double_vector &desc) = 0;
    virtual ~DescriptorKind();
    /*!
     */
     // ***************************************************************************************************************
     // Specialized descriptor initializing overloads
     // **************************************************************************************************************

    // ********Symmetry Functions********
    static DescriptorKind *
    initDescriptor(AvailableDescriptor availableDescriptorKind, std::vector<std::string> *species,
                   std::string *cutoff_function, double *cutoff_matrix,
                   std::vector<std::string> *symmetry_function_types, std::vector<int> *symmetry_function_sizes,
                   std::vector<double> *symmetry_function_parameters);

    // ********Bispectrum********
//    static DescriptorKind *
//    initDescriptor(AvailableDescriptor availableDescriptorKind, double rfac0_in, int twojmax_in, int diagonalstyle_in,
//                   int use_shared_arrays_in, double rmin0_in, int switch_flag_in, int bzero_flag_in,
//                   double * cutoff_array, std::vector<std::string> * species, std::vector<double> * weights);

};


#endif // __DESCRIPTOR_HPP_