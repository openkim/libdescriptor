#include "Descriptors.hpp"
#include "SymmetryFunctions/SymmetryFunctions.hpp"
#include "Bispectrum/Bispectrum.hpp"
#include "finite_difference.hpp"
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>

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
//void __enzyme_fwddiff(void (*) (int, int *, int *, int *, double *, double *, DescriptorKind *),
//                       int, int /* n_atoms */,
//                       int, int * /* Z */,
//                       int, int * /* neighbor list */,
//                       int, int * /* number_of_neigh_list */,
//                       int, double * /* coordinates */, double * /* derivative w.r.t coordinates */,
//                       int, double * /* zeta */, double * /* dzeta_dE */,
//                       int, DescriptorKind * /* DescriptorKind to diff */);

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
        default:
            std::cerr << "Descriptor kind not supported\n";
            throw std::invalid_argument("Descriptor kind not supported");
    }
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
    for(int i = 0; i < desc_kind->width; i++){dx_ddesc[i] = 0;}

    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 0, n_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++){
        d_coordinates[0] += dE_ddesc[i] * dx_ddesc[i];
    }
    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 1, n_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++){
        d_coordinates[1] += dE_ddesc[i] * dx_ddesc[i];
    }
    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 2, n_atoms * 3, desc_kind->width, dx_ddesc);
    for (int i = 0; i < desc_kind->width; i++){
        d_coordinates[2] += dE_ddesc[i] * dx_ddesc[i];
    }

    delete[] dx_ddesc;
}


DescriptorKind::~DescriptorKind() = default;