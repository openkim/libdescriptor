#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <fstream>
#include <vector>
#include "Descriptors.hpp"

#define COMMA ,

namespace py = pybind11;
using namespace Descriptor;

    class PyDescriptor: public Descriptor::DescriptorKind {
    public:
        using Descriptor::DescriptorKind::DescriptorKind;

            void compute(int index,
                         int n_atoms,
                         int *species,
                         int *neighbor_lists,
                         int number_of_neighbors,
                         double *coordinates,
                         double *desc) override {
            PYBIND11_OVERRIDE_PURE( void , DescriptorKind , compute,
                    index, n_atoms, species, neighbor_lists, number_of_neighbors, coordinates, desc
                    );
        };
    };

PYBIND11_MODULE(libdescriptor, m) {
    m.doc() = "Python interface to libdescriptor";

    py::enum_<AvailableDescriptor>(m, "AvailableDescriptors")
            .value("Bispectrum", AvailableDescriptor::KindBispectrum)
            .value("SymmetryFunctions", AvailableDescriptor::KindSymmetryFunctions);



//    py::class_<DescriptorKind>(m, "DescriptorKind")
//            .def(py::init<>());


}


