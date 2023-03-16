#include "Descriptors.hpp"
#include "SymmetryFunctions/SymmetryFunctions.hpp"
#include "Bispectrum/Bispectrum.hpp"
#include "finite_difference.hpp"
#include <vector>
#include <string>
#include <stdexcept>


using namespace Descriptor;

DescriptorKind *DescriptorKind::initDescriptor(std::string &descriptor_file_name,
                                               AvailableDescriptor descriptor_kind) {
    if (descriptor_kind == KindSymmetryFunctions) {
        auto sf = new SymmetryFunctions(descriptor_file_name);
        sf->descriptor_kind = descriptor_kind;
        sf->descriptor_param_file = descriptor_file_name;
        return sf;
    } else if (descriptor_kind == KindBispectrum) {
        auto bs = new Bispectrum(descriptor_file_name);
        bs->descriptor_kind = descriptor_kind;
        bs->descriptor_param_file = descriptor_file_name;
        return bs;
    } else {
        throw std::invalid_argument("Descriptor kind not implemented yet");
    }
}

DescriptorKind *DescriptorKind::initDescriptor(AvailableDescriptor descriptor_kind) {
    if (descriptor_kind == KindSymmetryFunctions) {
        return new SymmetryFunctions();
    } else if (descriptor_kind == KindBispectrum) {
        return new Bispectrum();
    } else {
        throw std::invalid_argument("Descriptor kind not implemented yet");
    }
}

/* This is wrapper for generic descriptor class compute function
   This wrapper will be differentiated using enzyme for more generic
   library structure. */
void Descriptor::compute(int const n_contributing_atoms /* contributing */,
                         int const n_total_atoms /* total */,
                         int *const species,
                         int *const neighbor_list,
                         int *const number_of_neighbors,
                         double *const coordinates,
                         double *const desc,
                         DescriptorKind *const desc_kind) {
    int *neighbor_ptr = neighbor_list;

    double_vector desc_vec(desc_kind->width * n_contributing_atoms);
    double_vector coordinates_vec(n_total_atoms * 3);

    for (int i = 0; i < n_total_atoms * 3; i++) {
        coordinates_vec(i) = coordinates[i];
    }
    double_vector desc_tmp(desc_kind->width);
    for (int i = 0; i < n_contributing_atoms; i++) {
        for (int j = 0; j < desc_kind->width; j++) {desc_tmp(j) = 0;} //init to zero
        desc_kind->compute(i, n_contributing_atoms, species, neighbor_ptr, number_of_neighbors[i],
                           coordinates_vec, desc_tmp);
        neighbor_ptr += number_of_neighbors[i];
        //copy to desc
        for (int j = 0; j < desc_kind->width; j++) {
            desc_vec(i * desc_kind->width + j) = desc_tmp(j);
        }
    }
    for (int i = 0; i < n_contributing_atoms * desc_kind->width; i++) {
        desc[i] = static_cast<double >(desc_vec(i));
    }
}

void Descriptor::gradient(int n_contributing_atoms /* contributing */,
                          int n_total_atoms /* total */,
                          int *species,
                          int *neighbor_list,
                          int *number_of_neighbors,
                          double *coordinates,
                          double *d_coordinates,
                          double *desc,
                          double *d_desc, /* vector for vjp or jvp */
                          DescriptorKind *desc_kind) {
    int *neighbor_ptr = neighbor_list;

//    double_vector desc_vec(desc_kind->width * n_contributing_atoms);
    double_vector coordinates_vec(n_total_atoms * 3);
//    double_vector d_coordinates_vec(n_total_atoms * 3);
    double_vector d_desc_vec(desc_kind->width * n_contributing_atoms);

    for (int i = 0; i < n_total_atoms * 3; i++) {
        coordinates_vec(i) = coordinates[i];
//        d_coordinates_vec(i) = d_coordinates[i];
    }

    for (int i = 0; i < desc_kind->width * n_contributing_atoms; i++) {
        d_desc_vec(i) = d_desc[i];
//        desc_vec(i) = desc[i];
    }

    auto wrapper_func = [&]( double_vector &coordinates_vec) {
        double_vector d_desc_vec_tmp(desc_kind->width * n_contributing_atoms);
        double_vector tmp_vec(desc_kind->width);
        for (int i = 0; i < n_contributing_atoms; i++) {
            for (int j = 0; j < desc_kind->width; j++) {tmp_vec(j) = 0;} //init to zero
            desc_kind->compute(i, n_contributing_atoms, species, neighbor_ptr, number_of_neighbors[i],
                               coordinates_vec, tmp_vec);
            neighbor_ptr += number_of_neighbors[i];
            //copy to desc
            for (int j = 0; j < desc_kind->width; j++) {
                d_desc_vec_tmp(i * desc_kind->width + j) = tmp_vec(j);
            }
        }
        return d_desc_vec_tmp;
    };

    double_matrix J = autodiff::jacobian(wrapper_func, wrt(coordinates_vec), at(coordinates_vec));
    double_vector d_coordinates_vec = d_desc_vec.transpose() * J;
    for (int i = 0; i < n_total_atoms * 3; i++) {
        d_coordinates[i] = static_cast<double>(d_coordinates_vec(i));
    }

}

void Descriptor::compute_single_atom(int index,
                                     int const n_total_atoms,
                                     int *const species,
                                     int *const neighbor_list,
                                     int number_of_neighbors,
                                     double *const coordinates,
                                     double *const desc,
                                     DescriptorKind *const desc_kind) {
    double_vector coordinates_vec(n_total_atoms * 3);
    for (int i = 0; i < n_total_atoms * 3; i++) {
        coordinates_vec(i) = coordinates[i];
    }
    double_vector desc_tmp(desc_kind->width);
    for (int j = 0; j < desc_kind->width; j++) {desc_tmp(j) = 0;} //init to zero
    desc_kind->compute(index, n_total_atoms, species, neighbor_list, number_of_neighbors,
                       coordinates_vec, desc_tmp);
    for (int j = 0; j < desc_kind->width; j++) {
        desc[j] = static_cast<double>(desc_tmp(j));
    }
}


void Descriptor::gradient_single_atom(int index,
                                      int n_total_atoms /* total */,
                                      int *species,
                                      int *neighbor_list,
                                      int number_of_neighbors,
                                      double *coordinates,
                                      double *d_coordinates,
                                      double *desc,
                                      double *d_desc, /* vector for vjp or jvp */
                                      DescriptorKind *desc_kind) {
    double_vector coordinates_vec(n_total_atoms * 3);
//    double_vector d_coordinates_vec(n_total_atoms * 3);
//    double_vector desc_vec(desc_kind->width);
    double_vector d_desc_vec(desc_kind->width);

    for (int i = 0; i < n_total_atoms * 3; i++) {
        coordinates_vec(i) = coordinates[i];
//        d_coordinates_vec(i) = d_coordinates[i];
    }
    for (int i = 0; i < desc_kind->width; i++) {
//        desc_vec(i) = desc[i];
        d_desc_vec(i) = d_desc[i];
    }

    auto fwd_func = [&](double_vector &coordinates_vec) {
        double_vector desc_tmp(desc_kind->width);
        desc_kind->compute(index, n_total_atoms, species, neighbor_list, number_of_neighbors,
                           coordinates_vec, desc_tmp);
        return desc_tmp;
    };

//    autodiff::gradient(fwd_func, wrt(coordinates_vec), at(d_desc_vec));
    double_matrix J = autodiff::jacobian(fwd_func, wrt(coordinates_vec), at(coordinates_vec));

    double_vector d_coord_tmp = d_desc_vec.transpose() * J;
    for (int i = 0; i < n_total_atoms * 3; i++) {
        d_coordinates[i] = static_cast<double>(d_coord_tmp(i));
    }

}

void Descriptor::num_gradient_single_atom(int index,
                                          int n_contributing_atoms /* contributing */,
                                          int n_total_atoms /* total */,
                                          int *species,
                                          int *neighbor_list,
                                          int number_of_neighbors,
                                          double *coordinates,
                                          double *d_coordinates,
                                          double *dE_ddesc, /* vector for vjp*/
                                          DescriptorKind *desc_kind) {
    auto f = [&](double *x, double *y) {
        double_vector coordinates_vec(n_total_atoms * 3);
        for (int i = 0; i < n_total_atoms * 3; i++) {
            coordinates_vec(i) = x[i];
        }
        double_vector desc_vec(desc_kind->width);
        for (int j = 0; j < desc_kind->width; j++) {desc_vec(j) = 0;} //init to zero
        desc_kind->compute(index, n_contributing_atoms, species, neighbor_list, number_of_neighbors, coordinates_vec, desc_vec);
        for (int j = 0; j < desc_kind->width; j++) {
            y[j] = static_cast<double>(desc_vec(j));
        }
    };

    auto dx_ddesc = new double[desc_kind->width];
    for(int i = 0; i < desc_kind->width; i++){dx_ddesc[i] = 0;}

    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 0, n_total_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++){
        d_coordinates[0] += dE_ddesc[i] * dx_ddesc[i];
    }
    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 1, n_total_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++){
        d_coordinates[1] += dE_ddesc[i] * dx_ddesc[i];
    }
    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 2, n_total_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++){
        d_coordinates[2] += dE_ddesc[i] * dx_ddesc[i];
    }

    delete[] dx_ddesc;
}

DescriptorKind::~DescriptorKind() = default;

DescriptorKind *
DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, std::vector<std::string> *species,
                               std::string *cutoff_function, double *cutoff_matrix,
                               std::vector<std::string> *symmetry_function_types,
                               std::vector<int> *symmetry_function_sizes,
                               std::vector<double> *symmetry_function_parameters){
    auto return_pointer = new SymmetryFunctions(species, cutoff_function, cutoff_matrix,
                                                   symmetry_function_types, symmetry_function_sizes,
                                                   symmetry_function_parameters);
    return return_pointer;
}

DescriptorKind *
DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, double rfac0_in, int twojmax_in, int diagonalstyle_in,
                   int use_shared_arrays_in, double rmin0_in, int switch_flag_in, int bzero_flag_in,
                   double * cutoff_array, std::vector<std::string> * species, std::vector<double> * weights){
    auto return_pointer = new Bispectrum(rfac0_in, twojmax_in, diagonalstyle_in,
                                            use_shared_arrays_in, rmin0_in, switch_flag_in, bzero_flag_in);
    return_pointer->width = return_pointer->get_width();
    return_pointer->set_species(species->size());
    std::string cutoff_function = "cos";
    return_pointer->set_cutoff(cutoff_function.c_str(), species->size(), cutoff_array);
    return_pointer->set_weight(species->size(), weights->data());
    return_pointer->descriptor_kind = availableDescriptorKind;
    return return_pointer;
}
