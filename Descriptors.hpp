// Generic header for inclusion and library design

#ifndef __DESCRIPTOR_HPP_
#define __DESCRIPTOR_HPP_

#include <map>
#include <string>
#include <vector>


namespace Descriptor {
    enum AvailableDescriptor {
        KindSymmetryFunctions, KindBispectrum
    };

    class DescriptorKind;
// TODO
//    void fwd_diff(int /* n_atoms */, int * /* Z */, int * /* neighbor list */, int * /* number_of_neigh_list */,
//                  double * /* coordinates */, double * /* d_coordinates */, double * /* zeta */,
//                  double * /* dE_dzeta */, DescriptorKind * /* DescriptorKind to diff */);

    void rev_diff(int /* n_atoms */, int * /* Z */, int * /* neighbor list */, int * /* number_of_neigh_list */,
                  double * /* coordinates */, double * /* d_coordinates */, double * /* zeta */,
                  double * /* dE_dzeta */, DescriptorKind * /* DescriptorKind to diff */);
// TODO
//    void jacobian(int /* n_atoms */, int * /* Z */, int * /* neighbor list */, int * /* number_of_neigh_list */,
//                  double * /* coordinates */, double * /* d_coordinates */, double * /* zeta */,
//                  double * /* dzeta_dr */, DescriptorKind * /* DescriptorKind to diff */);

    void compute(int /* n_atoms */, int * /* Z */, int * /* neighbor list */, int * /* number_of_neigh_list */,
                 double * /* coordinates */, double * /* zeta */,
                 DescriptorKind * /* DescriptorKind to diff */);
}


class Descriptor::DescriptorKind {
    // Base class for all descriptors
    // This will be the parent class for all descriptors. To ensure compatibility with
    // Enzyme, class members and datastructures will be kept simple. Structures like lists
    // of lists are excruciatingly slow to compile through using enzyme. Enzyme plays well with
    // more C-like code
public:
    AvailableDescriptor descriptor_kind;
    std::string descriptor_param_file;
    int width;

    DescriptorKind() = default;

    static DescriptorKind *initDescriptor(AvailableDescriptor);

    static DescriptorKind *initDescriptor(std::string &/* Filename */,
                                          AvailableDescriptor /* Descriptor kind */);

    virtual void compute(int /* index */,
                         int /* n_atoms */,
                         int * /* Z */,
                         int * /* neighbor lists */,
                         int  /* number of neighbors */,
                         double * /* coordinates */,
                         double * /* zeta */) = 0;

    ~DescriptorKind();
};


#endif // __DESCRIPTOR_HPP_