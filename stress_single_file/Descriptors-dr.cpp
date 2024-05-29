#include "Descriptors-dr.hpp"
//#include "Bispectrum.hpp"
//#include "SOAP.hpp"
//#include "Xi.hpp"
//#include "finite_difference.hpp"
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

using namespace Descriptor_dr;

DescriptorKind *DescriptorKind::initDescriptor(std::string &descriptor_file_name,
                                               AvailableDescriptor descriptor_kind) {
    if (descriptor_kind == KindSymmetryFunctions) {
        auto sf = new SymmetryFunctions_dr(descriptor_file_name);
        sf->descriptor_kind = descriptor_kind;
        sf->descriptor_param_file = descriptor_file_name;
        return sf;
//    } else if (descriptor_kind == KindBispectrum) {
//        auto bs = new Bispectrum(descriptor_file_name);
//        bs->descriptor_kind = descriptor_kind;
//        bs->descriptor_param_file = descriptor_file_name;
//        return bs;
//    } else if (descriptor_kind == KindSOAP) {
//        auto global = new SOAP(descriptor_file_name);
//        global->descriptor_kind = descriptor_kind;
//        global->descriptor_param_file = descriptor_file_name;
//        return global;
//    } else if (descriptor_kind == KindXi) {
//        auto xi = new Xi(descriptor_file_name);
//        xi->descriptor_kind = descriptor_kind;
//        xi->descriptor_param_file = descriptor_file_name;
//        return xi;
    } else {
        throw std::invalid_argument("Descriptor_dr kind not implemented yet");
    }
}

DescriptorKind *DescriptorKind::initDescriptor(AvailableDescriptor descriptor_kind) {
    if (descriptor_kind == KindSymmetryFunctions) {
        return new SymmetryFunctions_dr();
//    } else if (descriptor_kind == KindBispectrum) {
//        return new Bispectrum();
//    } else if (descriptor_kind == KindSOAP) {
//        return new SOAP();
//    } else if (descriptor_kind == KindXi) {
//        return new Xi();
    } else {
        throw std::invalid_argument("Descriptor_dr kind not implemented yet");
    }
}

/* This is wrapper for generic descriptor class compute function
   This wrapper will be differentiated using enzyme for more generic
   library structure. */
void Descriptor_dr::compute(int const n_atoms /* contributing */,
                            int *const species,
                            int *const neighbor_list,
                            int *const number_of_neighbors,
                            double *const distances,
                            double *const desc,
                            DescriptorKind *const desc_kind) {
    int *neighbor_ptr = neighbor_list;
    double *desc_ptr = desc;
    for (int i = 0; i < n_atoms; i++) {
        desc_kind->compute(i, n_atoms, species, neighbor_ptr, number_of_neighbors[i],
                           distances, desc_ptr);
        neighbor_ptr += number_of_neighbors[i];
        desc_ptr += desc_kind->width;
    }
}

void Descriptor_dr::gradient(int n_atoms /* contributing */,
                             int *species,
                             int *neighbor_list,
                             int *number_of_neighbors,
                             double *distances,
                             double *d_distances,
                             double *desc,
                             double *d_desc, /* vector for vjp or jvp */
                          DescriptorKind *desc_kind) {
    switch (desc_kind->descriptor_kind) {
        case KindSymmetryFunctions: {
            auto d_desc_kind = new SymmetryFunctions_dr();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));

            d_desc_kind->clone_empty(desc_kind);
            __enzyme_autodiff(compute, /* fn to be differentiated */
                              enzyme_const, n_atoms, /* Do not diff. against integer params */
                              enzyme_const, species,
                              enzyme_const, neighbor_list,
                              enzyme_const, number_of_neighbors,
                              enzyme_dup, distances, d_distances,
                              enzyme_dup, desc, d_desc,
                              enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
//        case KindBispectrum: {
//            auto d_desc_kind = new Bispectrum();
//            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
//
//            d_desc_kind->clone_empty(desc_kind);
//            __enzyme_autodiff(compute, /* fn to be differentiated */
//                              enzyme_const, n_atoms, /* Do not diff. against integer params */
//                              enzyme_const, species,
//                              enzyme_const, neighbor_list,
//                              enzyme_const, number_of_neighbors,
//                              enzyme_dup, coordinates, d_coordinates,
//                              enzyme_dup, desc, d_desc,
//                              enzyme_dup, desc_kind, d_desc_kind);
//            delete d_desc_kind;
//            return;
//        }
//        case KindSOAP: {
//            auto d_desc_kind = new SOAP();
//            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
//
//            d_desc_kind->clone_empty(desc_kind);
//            __enzyme_autodiff(compute, /* fn to be differentiated */
//                              enzyme_const, n_atoms, /* Do not diff. against integer params */
//                              enzyme_const, species,
//                              enzyme_const, neighbor_list,
//                              enzyme_const, number_of_neighbors,
//                              enzyme_dup, coordinates, d_coordinates,
//                              enzyme_dup, desc, d_desc,
//                              enzyme_dup, desc_kind, d_desc_kind);
//            // int *neighbor_ptr = neighbor_list;
//            // double *desc_ptr = desc;
//            // double *d_desc_ptr = d_desc;
//            // for (int i = 0; i < n_atoms; i++) {
//            //     __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
//            //                enzyme_const, i,
//            //                enzyme_const, n_atoms, /* Do not diff. against integer params */
//            //                enzyme_const, species,
//            //                enzyme_const, neighbor_ptr,
//            //                enzyme_const, number_of_neighbors[i],
//            //                enzyme_dup, coordinates, d_coordinates,
//            //                enzyme_dup, desc_ptr, d_desc_ptr,
//            //                enzyme_dup, desc_kind, d_desc_kind);
//            //     neighbor_ptr += number_of_neighbors[i];
//            //     desc_ptr += desc_kind->width;
//            //     d_desc_ptr += d_desc_kind->width;
//            // }
//            delete d_desc_kind;
//            return;
//        }
//        case KindXi: {
//            auto d_desc_kind = new Xi();
//            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
//
//            d_desc_kind->clone_empty(desc_kind);
//            __enzyme_autodiff(compute, /* fn to be differentiated */
//                              enzyme_const, n_atoms, /* Do not diff. against integer params */
//                              enzyme_const, species,
//                              enzyme_const, neighbor_list,
//                              enzyme_const, number_of_neighbors,
//                              enzyme_dup, coordinates, d_coordinates,
//                              enzyme_dup, desc, d_desc,
//                              enzyme_dup, desc_kind, d_desc_kind);
//            delete d_desc_kind;
//            return;
//        }
        default:
            std::cerr << "Descriptor_dr kind not supported\n";
            throw std::invalid_argument("Descriptor_dr kind not supported");
    }
}

void Descriptor_dr::compute_single_atom(int index,
                                        int const n_atoms /* contributing */,
                                        int *const species,
                                        int *const neighbor_list,
                                        int number_of_neighbors,
                                        double *const distances,
                                        double *const desc,
                                        DescriptorKind *const desc_kind) {
    desc_kind->compute(index, n_atoms, species, neighbor_list, number_of_neighbors,
                       distances, desc);
}


void Descriptor_dr::gradient_single_atom(int index,
                                         int n_atoms /* contributing */,
                                         int *species,
                                         int *neighbor_list,
                                         int number_of_neighbors,
                                         double *distances,
                                         double *d_distances,
                                         double *desc,
                                         double *d_desc, /* vector for vjp or jvp */
                                      DescriptorKind *desc_kind) {
    switch (desc_kind->descriptor_kind) {
        case KindSymmetryFunctions: {
            auto d_desc_kind = new SymmetryFunctions_dr();
            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
            d_desc_kind->clone_empty(desc_kind);

            __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
                                       enzyme_const, index,
                                       enzyme_const, n_atoms, /* Do not diff. against integer params */
                                       enzyme_const, species,
                                       enzyme_const, neighbor_list,
                                       enzyme_const, number_of_neighbors,
                                       enzyme_dup, distances, d_distances,
                                       enzyme_dup, desc, d_desc,
                                       enzyme_dup, desc_kind, d_desc_kind);
            delete d_desc_kind;
            return;
        }
//        case KindBispectrum: {
//            auto d_desc_kind = new Bispectrum();
//
//            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
//            d_desc_kind->clone_empty(desc_kind);
//
//            __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
//                                       enzyme_const, index,
//                                       enzyme_const, n_atoms, /* Do not diff. against integer params */
//                                       enzyme_const, species,
//                                       enzyme_const, neighbor_list,
//                                       enzyme_const, number_of_neighbors,
//                                       enzyme_dup, coordinates, d_coordinates,
//                                       enzyme_dup, desc, d_desc,
//                                       enzyme_dup, desc_kind, d_desc_kind);
//            delete d_desc_kind;
//            return;
//        }
//        case KindSOAP: {
//            auto d_desc_kind = new SOAP();
//
//            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
//            d_desc_kind->clone_empty(desc_kind);
//            __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
//                                       enzyme_const, index,
//                                       enzyme_const, n_atoms, /* Do not diff. against integer params */
//                                       enzyme_const, species,
//                                       enzyme_const, neighbor_list,
//                                       enzyme_const, number_of_neighbors,
//                                       enzyme_dup, coordinates, d_coordinates,
//                                       enzyme_dup, desc, d_desc,
//                                       enzyme_dup, desc_kind, d_desc_kind);
//            delete d_desc_kind;
//            return;
//        }
//        case KindXi: {
//            auto d_desc_kind = new Xi();
//
//            *((void **) d_desc_kind) = __enzyme_virtualreverse(*((void **) d_desc_kind));
//            d_desc_kind->clone_empty(desc_kind);
//            __enzyme_autodiff_one_atom(compute_single_atom, /* fn to be differentiated */
//                                       enzyme_const, index,
//                                       enzyme_const, n_atoms, /* Do not diff. against integer params */
//                                       enzyme_const, species,
//                                       enzyme_const, neighbor_list,
//                                       enzyme_const, number_of_neighbors,
//                                       enzyme_dup, coordinates, d_coordinates,
//                                       enzyme_dup, desc, d_desc,
//                                       enzyme_dup, desc_kind, d_desc_kind);
//            delete d_desc_kind;
//            return;
//        }
        default:
            std::cerr << "Descriptor_dr kind not supported\n";
            throw std::invalid_argument("Descriptor_dr kind not supported");
    }
}

void Descriptor_dr::jacobian(int n_atoms, /* contributing */
                          int n_total_atoms,
                             int *species,
                             int *neighbor_list,
                             int *number_of_neighs,
                             double *coordinates,
                             double *J_coordinates,
                             DescriptorKind *descriptor_to_diff) {
//    auto desc = new double[descriptor_to_diff->width * n_atoms];
//    std::vector<double> d_desc;
//    d_desc = std::vector<double>(descriptor_to_diff->width * n_atoms);
//    std::fill(desc, desc + descriptor_to_diff->width * n_atoms, 0.0);
//    // TODO: #pragma omp parallel for private(d_desc, desc)
//    for (int j = 0; j < descriptor_to_diff->width * n_atoms; j++) {
//        if ( j == 0 ) {
//            d_desc[j] = 1.0;
//        } else {
//            d_desc[j - 1] = 0.0;
//            d_desc[j] = 1.0;
//        }
//        gradient(n_atoms,
//                 species,
//                 neighbor_list,
//                 number_of_neighs,
//                 coordinates,
//                 J_coordinates + j * 3 * n_total_atoms,
//                 desc,
//                 d_desc.data(),
//                 descriptor_to_diff);
//    }
//
//    delete[] desc;
std::cerr << "Jacobian not implemented yet\n";
}

//void Descriptor_dr::num_gradient_single_atom(int index,
//                                             int n_atoms /* contributing */,
//                                             int *species,
//                                             int *neighbor_list,
//                                             int number_of_neighbors,
//                                             double *coordinates,
//                                             double *d_coordinates,
//                                             double *dE_ddesc, /* vector for vjp*/
//                                          DescriptorKind *desc_kind) {
//    auto f = [&](double *x, double *y) {
//        desc_kind->compute(index, n_atoms, species, neighbor_list, number_of_neighbors, x, y);
//    };
//
//    auto dx_ddesc = new double[desc_kind->width];
//    for (int i = 0; i < desc_kind->width; i++) { dx_ddesc[i] = 0; }
//
//    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 0, n_atoms * 3, desc_kind->width, dx_ddesc);
//    for (int i = 0; i < desc_kind->width; i++) {
//        d_coordinates[0] += dE_ddesc[i] * dx_ddesc[i];
//    }
//    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 1, n_atoms * 3, desc_kind->width, dx_ddesc);
//    for (int i = 0; i < desc_kind->width; i++) {
//        d_coordinates[1] += dE_ddesc[i] * dx_ddesc[i];
//    }
//    numdiff::vec_finite_difference_derivative(f, coordinates, index * 3 + 2, n_atoms * 3, desc_kind->width, dx_ddesc);
//    for (int i = 0; i < desc_kind->width; i++) {
//        d_coordinates[2] += dE_ddesc[i] * dx_ddesc[i];
//    }
//
//    delete[] dx_ddesc;
//}

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
//            return_pointer = new SymmetryFunctions_dr(species, cutoff_function, cutoff_matrix,
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
//            std::cerr << "Descriptor_dr kind not supported\n";
//            throw std::invalid_argument("Descriptor_dr kind not supported");
//    }
//
//    va_end(args);
//    return return_pointer;
//}

DescriptorKind *
DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, std::vector<std::string> *species,
                               std::string *cutoff_function, double *cutoff_matrix,
                               std::vector<std::string> *symmetry_function_types,
                               std::vector<int> *symmetry_function_sizes,
                               std::vector<double> *symmetry_function_parameters) {
    auto return_pointer = new SymmetryFunctions_dr(species, cutoff_function, cutoff_matrix,
                                                symmetry_function_types, symmetry_function_sizes,
                                                symmetry_function_parameters);
    return return_pointer;
}

//DescriptorKind *
//DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, double rfac0_in, int twojmax_in,
//                               int diagonalstyle_in,
//                               int use_shared_arrays_in, double rmin0_in, int switch_flag_in, int bzero_flag_in,
//                               double *cutoff_array, std::vector<std::string> *species, std::vector<double> *weights) {
//    auto return_pointer = new Bispectrum(rfac0_in, twojmax_in, diagonalstyle_in,
//                                         use_shared_arrays_in, rmin0_in, switch_flag_in, bzero_flag_in);
//    return_pointer->width = return_pointer->get_width();
//    return_pointer->set_species(species->size());
//    std::string cutoff_function = "cos";
//    return_pointer->set_cutoff(cutoff_function.c_str(), species->size(), cutoff_array);
//    return_pointer->set_weight(species->size(), weights->data());
//    return_pointer->descriptor_kind = availableDescriptorKind;
//    return return_pointer;
//}
//
//DescriptorKind *
//DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, int n_max, int l_max, double cutoff,
//                               std::vector<std::string> &species, std::string radial_basis, double eta) {
//    auto return_pointer = new SOAP(n_max, l_max, cutoff, species, radial_basis, eta);
//    return_pointer->width = return_pointer->get_width();
//    return_pointer->descriptor_kind = availableDescriptorKind;
//    return return_pointer;
//}
//
//DescriptorKind *
//DescriptorKind::initDescriptor(AvailableDescriptor availableDescriptorKind, int q, double cutoff,
//                               std::vector<std::string> &species, std::string& radial_basis) {
//    auto return_pointer = new Xi(q, cutoff, species, radial_basis);
//    return_pointer->width = return_pointer->get_width();
//    return_pointer->descriptor_kind = availableDescriptorKind;
//    return return_pointer;
//}

#ifdef DIM
#undef DIM
#endif
#define DIM 3

#ifdef MY_PI
#undef MY_PI
#endif
#define MY_PI 3.1415926535897932

inline void SymmetryFunctions_dr::set_species(std::vector<std::string> &species) {
    species_.resize(species.size());
    std::copy(species.begin(), species.end(), species_.begin());
}

inline void SymmetryFunctions_dr::get_species(std::vector<std::string> &species) {
    species.resize(species_.size());
    std::copy(species_.begin(), species_.end(), species.begin());
}

inline int SymmetryFunctions_dr::get_num_species() { return species_.size(); }

void SymmetryFunctions_dr::set_cutoff(char const *name,
                                   std::size_t const Nspecies,
                                   double const *rcut_2D) {
    (void) name;   // to avoid unused warning
    rcut_2D_.resize(Nspecies, Nspecies, rcut_2D);
}

inline double SymmetryFunctions_dr::get_cutoff(int const iCode, int const jCode) {
    return rcut_2D_(iCode, jCode);
}

void SymmetryFunctions_dr::add_descriptor(char const *name,
                                       double const *values,
                                       int const row,
                                       int const col) {
    if (strcmp(name, "g1") == 0) { name_.push_back(1); }
    if (strcmp(name, "g2") == 0) { name_.push_back(2); }
    if (strcmp(name, "g3") == 0) { name_.push_back(3); }
    if (strcmp(name, "g4") == 0) { name_.push_back(4); }
    if (strcmp(name, "g5") == 0) { name_.push_back(5); }

    Array2D<double> params(row, col, values);
    params_.push_back(std::move(params));

    auto sum = std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
    starting_index_.push_back(sum);

    num_param_sets_.push_back(row);
    num_params_.push_back(col);

    if (strcmp(name, "g4") == 0 || strcmp(name, "g5") == 0) {
        has_three_body_ = true;
    }
}

inline double cut_cos(double const r, double const rcut) {
    return (r < rcut) ? 0.5 * (std::cos(MY_PI * r / rcut) + 1.0) : 0.0;
}

int SymmetryFunctions_dr::get_num_descriptors() {
    return std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
}

void sym_g5(double const zeta, double const lambda, double const eta, double const *r,
            double const *rcut, double &phi) {
    double const rij = r[0];
    double const rik = r[1];
    double const rcutij = rcut[0];
    double const rcutik = rcut[1];

    if (rij > rcutij || rik > rcutik) { phi = 0.0; }
    else {
        double const rjk = r[2];
        double const rijsq = rij * rij;
        double const riksq = rik * rik;
        double const rjksq = rjk * rjk;

        // index is the apex atom
        double const cos_ijk = (rijsq + riksq - rjksq) / (2 * rij * rik);
        double const base = 1.0 + lambda * cos_ijk;

        // prevent numerical instability (when lambda=-1 and cos_ijk=1)
        double const costerm = (base <= 0) ? 0.0 : std::pow(base, zeta);
        double const eterm = std::exp(-eta * (rijsq + riksq));
        phi = std::pow(2, 1 - zeta) * costerm * eterm * cut_cos(rij, rcutij)
              * cut_cos(rik, rcutik);
    }
}

void sym_g4(double const zeta, double const lambda, double const eta, double const *r,
            double const *rcut, double &phi) {
    double const rij = r[0];
    double const rik = r[1];
    double const rjk = r[2];
    double const rcutij = rcut[0];
    double const rcutik = rcut[1];
    double const rcutjk = rcut[2];

    if (rij > rcutij || rik > rcutik || rjk > rcutjk) { phi = 0.0; }
    else {
        double const rijsq = rij * rij;
        double const riksq = rik * rik;
        double const rjksq = rjk * rjk;

        // index is the apex atom
        double const cos_ijk = (rijsq + riksq - rjksq) / (2 * rij * rik);
        double const base = 1 + lambda * cos_ijk;

        // prevent numerical instability (when lambda=-1 and cos_ijk=1)
        double const costerm = (base <= 0) ? 0.0 : std::pow(base, zeta);
        double const eterm = std::exp(-eta * (rijsq + riksq + rjksq));

        phi = std::pow(2, 1 - zeta) * costerm * eterm * cut_cos(rij, rcutij)
              * cut_cos(rik, rcutik) * cut_cos(rjk, rcutjk);
    }
}

void sym_g3(double const kappa, double const r, double const rcut, double &phi) {
    phi = std::cos(kappa * r) * cut_cos(r, rcut);
}

void sym_g2(double const eta, double const Rs, double const r, double const rcut, double &phi) {
    phi = std::exp(-eta * (r - Rs) * (r - Rs)) * cut_cos(r, rcut);
}

void sym_g1(double const r, double const rcut, double &phi) {
    phi = cut_cos(r, rcut);
}

void SymmetryFunctions_dr::compute(int const index,
                                int const n_atoms,
                                int *const species,
                                int *const neigh_list,
                                int const number_of_neigh,
                                double *const distances,
                                double *const desc) {
    // prepare data
    int const iSpecies = species[index];
    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < number_of_neigh; ++jj) {
        // adjust index of particle neighbor
        int const j = neigh_list[jj];
        int const jSpecies = species[j];
        // cutoff between ij
        double rcutij = rcut_2D_(iSpecies, jSpecies);;
        double const rijmag = distances[index * n_atoms + j];

        // if particles index and j not interact
        if (rijmag > rcutij) { continue; }

        // Loop over descriptors
        // two-body descriptors
        for (std::size_t p = 0; p < name_.size(); ++p) {
            if (name_[p] != 1 && name_[p] != 2 && name_[p] != 3) {
                continue;
            }

            int idx = starting_index_[p];
            // Loop over same descriptor but different parameter set
            for (int q = 0; q < num_param_sets_[p]; ++q) {
                double gc = 0.0;

                if (name_[p] == 1) {
                    sym_g1(rijmag, rcutij, gc);
                } else if (name_[p] == 2) {
                    double eta = params_[p](q, 0);
                    auto Rs = params_[p](q, 1);
                    sym_g2(eta, Rs, rijmag, rcutij, gc);
                } else if (name_[p] == 3) {
                    double kappa = params_[p](q, 0);
                    sym_g3(kappa, rijmag, rcutij, gc);
                }
                desc[idx] += gc;
                ++idx;
            }
        }
        // three-body descriptors
        if (has_three_body_ == 0) { continue; }

        // Loop over kk
        for (int kk = jj + 1; kk < number_of_neigh; ++kk) {
            // Adjust index of particle neighbor
            int const k = neigh_list[kk];
            int const kSpecies = species[k];

            // cutoff between ik and jk
            double const rcutik = rcut_2D_[iSpecies][kSpecies];
            double const rcutjk = rcut_2D_[jSpecies][kSpecies];

            // Compute rik, rjk and their squares
            double const rikmag = distances[index * n_atoms + k];
            double const rjkmag = distances[j * n_atoms + k];

            // Check whether three-dody not interacting
            if (rikmag > rcutik) { continue; }

            double const rvec[3] = {rijmag, rikmag, rjkmag};
            double const rcutvec[3] = {rcutij, rcutik, rcutjk};

            // Loop over descriptors
            // three-body descriptors
            for (size_t p = 0; p < name_.size(); ++p) {
                if (name_[p] != 4 && name_[p] != 5) { continue; }
                int idx = starting_index_[p];

                // Loop over same descriptor but different parameter set
                for (int q = 0; q < num_param_sets_[p]; ++q) {
                    double gc = 0.0;

                    if (name_[p] == 4) {
                        double zeta = params_[p](q, 0);
                        double lambda = params_[p](q, 1);
                        double eta = params_[p](q, 2);

                        sym_g4(zeta, lambda, eta, rvec, rcutvec, gc);
                    } else if (name_[p] == 5) {
                        double zeta = params_[p](q, 0);
                        double lambda = params_[p](q, 1);
                        double eta = params_[p](q, 2);

                        sym_g5(zeta, lambda, eta, rvec, rcutvec, gc);
                    }

                    desc[idx] += gc;
                    ++idx;
                }
            }
        }
    }
}


SymmetryFunctions_dr::SymmetryFunctions_dr(std::string &file_name) {
    initFromFile(file_name);
}

void SymmetryFunctions_dr::initFromFile(std::string &file_name) {

    // Open the descriptor file
    std::ifstream file = FileIOUtils::open_file(file_name);

    // String containing data line and list of parameters
    std::vector<std::string> string_params;
    std::vector<double> double_params;
    std::vector<int> int_params;
    std::vector<bool> bool_params;
    std::string line;

    // Params
    std::string cutoff_type;
    int num_species;
    std::vector<std::string> species;
    std::map<std::pair<std::string, std::string>, double> species_pair_cutoffs;
    bool normalize = false;
    int descriptor_size;
    std::vector<std::string> symmetry_function_types;
    std::map<std::string, std::vector<double>> symmetry_function_params;

    // read Cutoff type name
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_string_params(line, string_params, 1);
    cutoff_type = string_params[0];
    string_params.clear();
    line.clear();

    // read number of species
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_int_params(line, int_params, 1);
    num_species = int_params[0];
    int_params.clear();
    line.clear();

    // read species names and cutoffs
    std::set<std::string> species_set;
    std::pair<std::string, std::string> species_species_pair;
    int species_index = num_species;
    do {
        FileIOUtils::get_next_data_line(file, line);
        FileIOUtils::parse_string_params(line, string_params, 2);
        FileIOUtils::parse_double_params(line, double_params, 1);

        species_set.insert(string_params[0]);
        species_set.insert(string_params[1]);

        species_species_pair = std::make_pair(string_params[0], string_params[1]);
        species_pair_cutoffs[species_species_pair] = double_params[0];

        species_species_pair = std::make_pair(string_params[1], string_params[0]);
        species_pair_cutoffs[species_species_pair] = double_params[0];

        double_params.clear();
        string_params.clear();
        line.clear();
    } while (--species_index > 0);
    species = std::vector<std::string>(species_set.begin(), species_set.end());

    // number of symmetry function types
    int num_symmetry_function_types;
    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_int_params(line, int_params, 1);
    num_symmetry_function_types = int_params[0];
    int_params.clear();
    line.clear();

    // Read symmetry function types and their parameters
    std::vector<double> symmetry_function_params_vector;
    std::map<std::string, std::pair<int, int>> symmetry_function_type_param_sizes;

    for (int i = 0; i < num_symmetry_function_types; i++) {
        FileIOUtils::get_next_data_line(file, line);
        FileIOUtils::parse_string_params(line, string_params, 1);
        symmetry_function_types.push_back(string_params[0]);
        string_params.clear();

        FileIOUtils::parse_int_params(line, int_params, 2);
        symmetry_function_type_param_sizes[symmetry_function_types[i]] = std::make_pair(int_params[0], int_params[1]);
        int_params.clear();
        line.clear();

        for (int j = 0; j < symmetry_function_type_param_sizes[symmetry_function_types[i]].first; j++) {
            FileIOUtils::get_next_data_line(file, line);
            FileIOUtils::parse_double_params(line, double_params,
                                             symmetry_function_type_param_sizes[symmetry_function_types[i]].second);
            for (int k = 0; k < symmetry_function_type_param_sizes[symmetry_function_types[i]].second; k++) {
                symmetry_function_params_vector.push_back(double_params[k]);
            }
            symmetry_function_params[symmetry_function_types[i]] = symmetry_function_params_vector;

            double_params.clear();
            line.clear();
        }
    }

    int_params.clear();
    double_params.clear();
    string_params.clear();
    line.clear();

    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_bool_params(line, bool_params, 1);
    normalize = bool_params[0];
    bool_params.clear();
    line.clear();

    FileIOUtils::get_next_data_line(file, line);
    FileIOUtils::parse_int_params(line, int_params, 1);
    descriptor_size = int_params[0];
    int_params.clear();
    line.clear();

    file.close();

    // Initialize the descriptor
    auto cutoff_matrix = std::make_unique<double[]>(num_species * num_species);
    for (int i = 0; i < num_species; i++) {
        for (int j = 0; j < num_species; j++) {
            species_species_pair = std::make_pair(species[i], species[j]);
            cutoff_matrix[i * num_species + j] = species_pair_cutoffs[species_species_pair];
        }
    }

    // set cut offs
    n_species = num_species;
    set_species(species);

    set_cutoff(cutoff_type.c_str(), num_species, cutoff_matrix.get());

    // set symmetry functions and their parameters and width
    int _width = 0;
    for (int i = 0; i < num_symmetry_function_types; i++) {
        add_descriptor(symmetry_function_types[i].c_str(),
                       symmetry_function_params[symmetry_function_types[i]].data(),
                       symmetry_function_type_param_sizes[symmetry_function_types[i]].first,
                       symmetry_function_type_param_sizes[symmetry_function_types[i]].second);
        _width += symmetry_function_type_param_sizes[symmetry_function_types[i]].first;
    }
    width = _width;
}


void SymmetryFunctions_dr::clone_empty(DescriptorKind *descriptorKind) {
    auto d_sf = dynamic_cast<SymmetryFunctions_dr *>(descriptorKind);
    name_ = d_sf->name_;
    params_ = d_sf->params_;
    rcut_2D_ = d_sf->rcut_2D_;
    has_three_body_ = d_sf->has_three_body_;
    width = d_sf->width;
    num_param_sets_ = d_sf->num_param_sets_;
    num_params_ = d_sf->num_params_;
    // set params to zero, to differentiate against
    for (int i = 0; i < name_.size(); i++) {
        for (int j = 0; j < num_param_sets_[i]; j++) {
            for (int k = 0; k < num_params_[i]; k++) {
                params_[i](j, k) = 0.0;
            }
        }
    }
}

SymmetryFunctions_dr::SymmetryFunctions_dr(std::vector<std::string> *species, std::string *cutoff_function,
                                     double *cutoff_matrix, std::vector<std::string> *symmetry_function_types,
                                     std::vector<int> *symmetry_function_sizes,
                                     std::vector<double> *symmetry_function_parameters) {
    // set cut offs
    n_species = species->size();
    set_species(*species);

    set_cutoff(cutoff_function->c_str(), n_species, cutoff_matrix);

    // set symmetry functions and their parameters
    int offset = 0;
    int _width = 0;
    for (int i = 0; i < symmetry_function_types->size(); i++) {
        add_descriptor(symmetry_function_types->at(i).c_str(),
                       symmetry_function_parameters->data() + offset,
                       symmetry_function_sizes->at(2 * i), symmetry_function_sizes->at(2 * i + 1));
        offset += symmetry_function_sizes->at(2 * i) * symmetry_function_sizes->at(2 * i + 1);
        _width += symmetry_function_sizes->at(2 * i);
    }
    width = _width;
}

