#include "Descriptors.hpp"
#include "SymmetryFunctions.hpp"
#include "Bispectrum.hpp"
#include "SOAP.hpp"
#include "Xi.hpp"
#include "finite_difference.hpp"
#include <vector>
#include <string>
#include <stdexcept>
//#include <omp.h>

int enzyme_dup, enzyme_out, enzyme_const;

template<typename T>
T __enzyme_virtualreverse(T);

// Rev mode diff
void __enzyme_autodiff(void (*)(int, int *, int *, int *, double *, double *, DescriptorKind *),
                       int, int /* n_atoms */,
                       int, int * /* Z */,
                       int, int * /* neighbor list */,
                       int, int * /* number_of_neigh_list */,
                       int, double * /* coordinates */, double * /* derivative w.r.t coordinates */,
                       int, double * /* zeta */, double * /* dzeta_dE */,
                       int, DescriptorKind * /* DescriptorKind to diff */, DescriptorKind * /* d_DescriptorKind */);

void __enzyme_autodiff_batch(void (*)(int, int *, int *, int *, int *, int *, double *, double *, DescriptorKind *),
                       int, int, /* n_configs */
                       int, int *, /* n_atoms */
                       int, int *, /* ptr */
                       int, int *, /* Z */
                       int, int *, /* neighbor list */
                       int, int *, /* number_of_neigh_list */
                       int, double * /* coordinates */, double * /* derivative w.r.t coordinates */,
                       int, double * /* zeta */, double * /* dzeta_dE */,
                       int, DescriptorKind * /* DescriptorKind to diff */, DescriptorKind * /* d_DescriptorKind */);

void __enzyme_autodiff_one_atom(void (*)(int, int, int *, int *, int, double *, double *, DescriptorKind *),
                                int, int,
                                int, int /* n_atoms */,
                                int, int * /* Z */,
                                int, int * /* neighbor list */,
                                int, int  /* number_of_neigh_list */,
                                int, double * /* coordinates */, double * /* derivative w.r.t coordinates */,
                                int, double * /* zeta */, double * /* dzeta_dE */,
                                int, DescriptorKind * /* DescriptorKind to diff */,
                                DescriptorKind * /* d_DescriptorKind */);

/* Fwd mode AD. As Descriptors are many-to-many functions, in certain situations,
 * like in case of global descriptors, FWD mode might be more performant.
 * Also in future when we will have full jacobian of descriptors enabled,
 * FWD mode might be used. So adding it now only as a simple to-be-used in future.
 */
//TODO
void __enzyme_fwddiff(void (*)(int, int *, int *, int *, double *, double *, DescriptorKind *),
                      int, int /* n_atoms */,
                      int, int * /* Z */,
                      int, int * /* neighbor list */,
                      int, int * /* number_of_neigh_list */,
                      int, double * /* coordinates */, double * /* derivative w.r.t coordinates */,
                      int, double * /* zeta */, double * /* dzeta_dE */,
                      int, DescriptorKind * /* DescriptorKind to diff */, DescriptorKind * /* d_DescriptorKind */);

void __enzyme_fwddiff_one_atom(void (*)(int, int, int *, int *, int, double *, double *, DescriptorKind *),
                               int, int,
                               int, int /* n_atoms */,
                               int, int * /* Z */,
                               int, int * /* neighbor list */,
                               int, int  /* number_of_neigh_list */,
                               int, double * /* coordinates */, double * /* derivative w.r.t coordinates */,
                               int, double * /* zeta */, double * /* dzeta_dE */,
                               int, DescriptorKind * /* DescriptorKind to diff */,
                               DescriptorKind * /* d_DescriptorKind */);

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
    } else if (descriptor_kind == KindSOAP) {
        auto global = new SOAP(descriptor_file_name);
        global->descriptor_kind = descriptor_kind;
        global->descriptor_param_file = descriptor_file_name;
        return global;
    } else if (descriptor_kind == KindXi) {
        auto xi = new Xi(descriptor_file_name);
        xi->descriptor_kind = descriptor_kind;
        xi->descriptor_param_file = descriptor_file_name;
        return xi;
    } else {
        throw std::invalid_argument("Descriptor kind not implemented yet");
    }
}

DescriptorKind *DescriptorKind::initDescriptor(AvailableDescriptor descriptor_kind) {
    if (descriptor_kind == KindSymmetryFunctions) {
        return new SymmetryFunctions();
    } else if (descriptor_kind == KindBispectrum) {
        return new Bispectrum();
    } else if (descriptor_kind == KindSOAP) {
        return new SOAP();
    } else if (descriptor_kind == KindXi) {
        return new Xi();
    } else {
        throw std::invalid_argument("Descriptor kind not implemented yet");
    }
}

/* This is wrapper for generic descriptor class compute function
   This wrapper will be differentiated using enzyme for more generic
   library structure. */
void Descriptor::compute(int const n_atoms /* contributing */,
                         int *const species,
                         int *const neighbor_list,
                         int *const number_of_neighbors,
                         double *const coordinates,
                         double *const desc,
                         DescriptorKind *const desc_kind) {
    int *neighbor_ptr = neighbor_list;
    double *desc_ptr = desc;
    for (int i = 0; i < n_atoms; i++) {
        desc_kind->compute(i, n_atoms, species, neighbor_ptr, number_of_neighbors[i],
                           coordinates, desc_ptr);
        neighbor_ptr += number_of_neighbors[i];
        desc_ptr += desc_kind->width;
    }
}

void Descriptor::compute_batch(int n_configurations, int *n_atoms, int *ptr, int *species, int *neighbor_list,
                               int *number_of_neighs, double *coordinates, double *desc,
                               Descriptor::DescriptorKind *descriptor_kind) {
    int *neighbor_ptr = neighbor_list;
    double *desc_ptr = desc;
    int idx = 0;
    int total_atoms = 0;
    for (int i = 0; i < n_configurations; i++) {
        total_atoms += n_atoms[i];
    }

    int total_neighs = 0;
    for (int i = 0; i < total_atoms; i++) {
        total_neighs += number_of_neighs[i];
    }

    for (int i = 0; i < n_configurations; i++) {
        for (int j = 0; j < n_atoms[i]; j++) {
            descriptor_kind->compute(j + ptr[i], total_atoms, species, neighbor_ptr, number_of_neighs[idx],
                                     coordinates, desc_ptr);
            neighbor_ptr += number_of_neighs[idx];
            desc_ptr += descriptor_kind->width;
            idx++;
        }
    }
}

void Descriptor::gradient(int n_atoms /* contributing */,
                          int *species,
                          int *neighbor_list,
                          int *number_of_neighbors,
                          double *coordinates,
                          double *d_coordinates,
                          double *desc,
                          double *d_desc, /* vector for vjp or jvp */
                          DescriptorKind *desc_kind) {
    switch (desc_kind->descriptor_kind) {
        case KindSymmetryFunctions: {
            auto d_desc_kind = new SymmetryFunctions();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(desc_kind);
            __enzyme_autodiff(compute, /* fn to be differentiated */
                              enzyme_const, n_atoms, /* Do not diff. against integer params */
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighbors,
                              enzyme_dup, coordinates, d_coordinates,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        case KindBispectrum: {
            auto d_desc_kind = new Bispectrum();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(desc_kind);
            __enzyme_autodiff(compute, /* fn to be differentiated */
                              enzyme_const, n_atoms, /* Do not diff. against integer params */
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighbors,
                              enzyme_dup, coordinates, d_coordinates,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        case KindSOAP: {
            auto d_desc_kind = new SOAP();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(desc_kind);
            __enzyme_autodiff(compute, /* fn to be differentiated */
                              enzyme_const, n_atoms, /* Do not diff. against integer params */
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighbors,
                              enzyme_dup, coordinates, d_coordinates,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, desc_kind, d_desc_kind);
            // int *neighbor_ptr = neighbor_list;
            // double *desc_ptr = desc;
            // double *d_desc_ptr = d_desc;
            // for (int i = 0; i < n_atoms; i++) {
            //     __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
            //                enzyme_const, i,
            //                enzyme_const, n_atoms, /* Do not diff. against integer params */
            //                enzyme_const, species,
            //                enzyme_const, neighbor_ptr,
            //                enzyme_const, number_of_neighbors[i],
            //                enzyme_dup, coordinates, d_coordinates,
            //                enzyme_dup, desc_ptr, d_desc_ptr,
            //                enzyme_dup, desc_kind, d_desc_kind);
            //     neighbor_ptr += number_of_neighbors[i];
            //     desc_ptr += desc_kind->width;
            //     d_desc_ptr += d_desc_kind->width;
            // }
            delete d_desc_kind;
            return;
        }
        case KindXi: {
            auto d_desc_kind = new Xi();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(desc_kind);
            __enzyme_autodiff(compute, /* fn to be differentiated */
                              enzyme_const, n_atoms, /* Do not diff. against integer params */
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighbors,
                              enzyme_dup, coordinates, d_coordinates,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        default:
            std::cerr << "Descriptor kind not supported\n";
            throw std::invalid_argument("Descriptor kind not supported");
    }
}


void Descriptor::gradient_batch(int n_configurations,
                                int *n_atoms,
                                int *ptr,
                                int *species,
                                int *neighbor_list,
                                int *number_of_neighs,
                                double *coordinates,
                                double *d_coordinates,
                                double *desc,
                                double *d_desc,
                                Descriptor::DescriptorKind *descriptor_to_diff) {
    switch (descriptor_to_diff->descriptor_kind) {
        case KindSymmetryFunctions: {
            auto d_desc_kind = new SymmetryFunctions();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(descriptor_to_diff);
            __enzyme_autodiff_batch(compute_batch, /* fn to be differentiated */
                              enzyme_const, n_configurations, /* Do not diff. against integer params */
                              enzyme_const, n_atoms,
                              enzyme_const, ptr,
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighs,
                              enzyme_dup, coordinates, d_coordinates,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, descriptor_to_diff, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        case KindBispectrum: {
            auto d_desc_kind = new Bispectrum();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(descriptor_to_diff);
            __enzyme_autodiff_batch(compute_batch, /* fn to be differentiated */
                              enzyme_const, n_configurations, /* Do not diff. against integer params */
                              enzyme_const, n_atoms,
                              enzyme_const, ptr,
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighs,
                              enzyme_dup, coordinates, d_coordinates,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, descriptor_to_diff, d_desc_kind);
            delete d_desc_kind;
            return;
            return;
        }
        case KindSOAP: {
            auto d_desc_kind = new SOAP();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(descriptor_to_diff);
            __enzyme_autodiff_batch(compute_batch, /* fn to be differentiated */
                              enzyme_const, n_configurations, /* Do not diff. against integer params */
                              enzyme_const, n_atoms,
                              enzyme_const, ptr,
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighs,
                              enzyme_dup, coordinates, d_coordinates,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, descriptor_to_diff, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        case KindXi: {
            auto d_desc_kind = new Xi();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(descriptor_to_diff);
            __enzyme_autodiff_batch(compute_batch, /* fn to be differentiated */
                              enzyme_const, n_configurations, /* Do not diff. against integer params */
                              enzyme_const, n_atoms,
                              enzyme_const, ptr,
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighs,
                              enzyme_dup, coordinates, d_coordinates,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, descriptor_to_diff, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        default:
            std::cerr << "Descriptor kind not supported\n";
            throw std::invalid_argument("Descriptor kind not supported");
    }
}


void Descriptor::compute_single_atom(int index,
                                     int const n_atoms /* contributing */,
                                     int *const species,
                                     int *const neighbor_list,
                                     int number_of_neighbors,
                                     double *const coordinates,
                                     double *const desc,
                                     DescriptorKind *const desc_kind) {
    desc_kind->compute(index, n_atoms, species, neighbor_list, number_of_neighbors,
                       coordinates, desc);
}


void Descriptor::gradient_single_atom(int index,
                                      int n_atoms /* contributing */,
                                      int *species,
                                      int *neighbor_list,
                                      int number_of_neighbors,
                                      double *coordinates,
                                      double *d_coordinates,
                                      double *desc,
                                      double *d_desc, /* vector for vjp or jvp */
                                      DescriptorKind *desc_kind) {
    switch (desc_kind->descriptor_kind) {
        case KindSymmetryFunctions: {
            auto d_desc_kind = new SymmetryFunctions();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
            d_desc_kind->clone_empty(desc_kind);

            __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
                                       enzyme_const, index,
                                       enzyme_const, n_atoms, /* Do not diff. against integer params */
                                       enzyme_const, species,
                                       enzyme_const, neighbor_list,
                                       enzyme_const, number_of_neighbors,
                                       enzyme_dup, coordinates, d_coordinates,
                                       enzyme_dup, desc, d_desc,
                                       enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        case KindBispectrum: {
            auto d_desc_kind = new Bispectrum();

            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
            d_desc_kind->clone_empty(desc_kind);

            __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
                                       enzyme_const, index,
                                       enzyme_const, n_atoms, /* Do not diff. against integer params */
                                       enzyme_const, species,
                                       enzyme_const, neighbor_list,
                                       enzyme_const, number_of_neighbors,
                                       enzyme_dup, coordinates, d_coordinates,
                                       enzyme_dup, desc, d_desc,
                                       enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        case KindSOAP: {
            auto d_desc_kind = new SOAP();

            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
            d_desc_kind->clone_empty(desc_kind);
            __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
                                       enzyme_const, index,
                                       enzyme_const, n_atoms, /* Do not diff. against integer params */
                                       enzyme_const, species,
                                       enzyme_const, neighbor_list,
                                       enzyme_const, number_of_neighbors,
                                       enzyme_dup, coordinates, d_coordinates,
                                       enzyme_dup, desc, d_desc,
                                       enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        case KindXi: {
            auto d_desc_kind = new Xi();

            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
            d_desc_kind->clone_empty(desc_kind);
            __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
                                       enzyme_const, index,
                                       enzyme_const, n_atoms, /* Do not diff. against integer params */
                                       enzyme_const, species,
                                       enzyme_const, neighbor_list,
                                       enzyme_const, number_of_neighbors,
                                       enzyme_dup, coordinates, d_coordinates,
                                       enzyme_dup, desc, d_desc,
                                       enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
        default:
            std::cerr << "Descriptor kind not supported\n";
            throw std::invalid_argument("Descriptor kind not supported");
    }
}

void Descriptor::jacobian(int n_atoms, /* contributing */
                          int n_total_atoms,
                          int *species,
                          int *neighbor_list,
                          int *number_of_neighs,
                          double *coordinates,
                          double *J_coordinates,
                          DescriptorKind *descriptor_to_diff) {
    auto desc = new double[descriptor_to_diff->width * n_atoms];
    std::vector<double> d_desc;
    d_desc = std::vector<double>(descriptor_to_diff->width * n_atoms);
    std::fill(desc, desc + descriptor_to_diff->width * n_atoms, 0.0);
    // TODO: #pragma omp parallel for private(d_desc, desc)
    for (int j = 0; j < descriptor_to_diff->width * n_atoms; j++) {
        if ( j == 0 ) {
            d_desc[j] = 1.0;
        } else {
            d_desc[j - 1] = 0.0;
            d_desc[j] = 1.0;
        }
        gradient(n_atoms,
                 species,
                 neighbor_list,
                 number_of_neighs,
                 coordinates,
                 J_coordinates + j * 3 * n_total_atoms,
                 desc,
                 d_desc.data(),
                 descriptor_to_diff);
    }

    delete[] desc;
}

void Descriptor::num_gradient_single_atom(int index,
                                          int n_atoms /* contributing */,
                                          int *species,
                                          int *neighbor_list,
                                          int number_of_neighbors,
                                          double *coordinates,
                                          double *d_coordinates,
                                          double *dE_ddesc, /* vector for vjp*/
                                          DescriptorKind *desc_kind) {
    auto f = [&](double *x, double *y) {
        desc_kind->compute(index, n_atoms, species, neighbor_list, number_of_neighbors, x, y);
    };

    auto dx_ddesc = new double[desc_kind->width];
    for (int i = 0; i < desc_kind->width; i++) { dx_ddesc[i] = 0; }

    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 0, n_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++) {
        d_coordinates[0] += dE_ddesc[i] * dx_ddesc[i];
    }
    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 1, n_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++) {
        d_coordinates[1] += dE_ddesc[i] * dx_ddesc[i];
    }
    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 2, n_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++) {
        d_coordinates[2] += dE_ddesc[i] * dx_ddesc[i];
    }

    delete[] dx_ddesc;
}

DescriptorKind::~DescriptorKind() = default;

// *********************************************************************************************************************
// Functions for individual descriptor kinds initialization, so that Pybind11 file can remain light and clean
// Revisit this design later with variadic templates if it is more feasible
// *********************************************************************************************************************

// TODO: URGENT: C++ compliant variadic ags:
//  https://wiki.sei.cmu.edu/confluence/display/cplusplus/DCL50-C++PP.+Do+not+define+a+C-style+variadic+function
// This needs to be done manually for now it seems. Pybind11 does not support variadic args.
// Therefore for Python bindings, every descriptor needs individual init function. This will be
// placed in a seperate file.
// TODO: Revisit in future with variadic template. Pybind11 might support it.
//DescriptorKind *DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, ...) {
//    va_list args;
//    va_start(args, availableDescriptorKind);
//
//    DescriptorKind * return_pointer = nullptr;
//
//    switch (availableDescriptorKind) {
//        case KindSymmetryFunctions: {
//            std::vector<std::string> *species;
//            std::string *cutoff_function;
//            double *cutoff_matrix;
//            std::vector<std::string> *symmetry_function_types;
//            std::vector<int> *symmetry_function_sizes;
//            std::vector<double> *symmetry_function_parameters;
//            species = va_arg(args, std::vector<std::string> *);
//            cutoff_function = va_arg(args, std::string *);
//            cutoff_matrix = va_arg(args, double *);
//            symmetry_function_types = va_arg(args, std::vector<std::string> *);
//            symmetry_function_sizes = va_arg(args, std::vector<int> *);
//            symmetry_function_parameters = va_arg(args, std::vector<double> *);
//
//            return_pointer = new SymmetryFunctions(species, cutoff_function, cutoff_matrix,
//                                                   symmetry_function_types, symmetry_function_sizes,
//                                                   symmetry_function_parameters);
//        }
//        case KindBispectrum: {
//            double rfac0_in = va_arg(args, double);
//            int twojmax_in = va_arg(args, int);
//            int diagonalstyle_in = va_arg(args, int);
//            int use_shared_arrays_in = va_arg(args, int);
//            double rmin0_in = va_arg(args, double);
//            int switch_flag_in = va_arg(args, int);
//            int bzero_flag_in = va_arg(args, int);
//
//            return_pointer = new Bispectrum(rfac0_in, twojmax_in, diagonalstyle_in,
//                                            use_shared_arrays_in, rmin0_in, switch_flag_in, bzero_flag_in);
//        }
//        default:
//            std::cerr << "Descriptor kind not supported\n";
//            throw std::invalid_argument("Descriptor kind not supported");
//    }
//
//    va_end(args);
//    return return_pointer;
//}

DescriptorKind *
DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, std::vector<std::string> *species,
                               std::string *cutoff_function, double *cutoff,
                               std::vector<std::string> *symmetry_function_types,
                               std::vector<int> *symmetry_function_sizes,
                               std::vector<double> *symmetry_function_parameters) {
    std::vector<double> cutoff_array(species->size() * species->size(), cutoff[0]);
    auto return_pointer = new SymmetryFunctions(species, cutoff_function, cutoff_array.data(),
                                                symmetry_function_types, symmetry_function_sizes,
                                                symmetry_function_parameters);
    return return_pointer;
}

DescriptorKind *
DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, double rfac0_in, int twojmax_in,
                               int diagonalstyle_in,
                               int use_shared_arrays_in, double rmin0_in, int switch_flag_in, int bzero_flag_in,
                               double *cutoff, std::vector<std::string> *species, std::vector<double> *weights) {
    std::vector<double> cutoff_array(species->size() * species->size(), cutoff[0]);
    auto return_pointer = new Bispectrum(rfac0_in, twojmax_in, diagonalstyle_in,
                                         use_shared_arrays_in, rmin0_in, switch_flag_in, bzero_flag_in);
    return_pointer->width = return_pointer->get_width();
    return_pointer->set_species(species->size());
    std::string cutoff_function = "cos";
    return_pointer->set_cutoff(cutoff_function.c_str(), species->size(), cutoff_array.data());
    return_pointer->set_weight(species->size(), weights->data());
    return_pointer ->descriptor_kind = availableDescriptorKind;
    return return_pointer;
}

DescriptorKind *
DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, int n_max, int l_max, double cutoff,
                               std::vector<std::string> &species, std::string radial_basis, double eta) {
    auto return_pointer = new SOAP(n_max, l_max, cutoff, species, radial_basis, eta);
    return_pointer->width = return_pointer->get_width();
    return_pointer->descriptor_kind = availableDescriptorKind;

    return return_pointer;
}

DescriptorKind *
DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, int q, double cutoff,
                               std::vector<std::string> &species, std::string& radial_basis) {
    auto return_pointer = new Xi(q, cutoff, species, radial_basis);
    return_pointer->width = return_pointer->get_width();
    return_pointer->descriptor_kind = availableDescriptorKind;
    return return_pointer;
}