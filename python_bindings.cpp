#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <fstream>
#include <vector>
#include <string>
#include "Descriptors.hpp"

namespace py = pybind11;
using namespace Descriptor;

class PyDescriptorKind : public DescriptorKind {
public:
    using DescriptorKind::DescriptorKind;

    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neighbor_lists,
                 int number_of_neighbors,
                 double *coordinates,
                 double *desc) override {
        PYBIND11_OVERRIDE_PURE(void, DescriptorKind, compute,
                               index, n_atoms, species, neighbor_lists, number_of_neighbors, coordinates, desc
        );
    };
};

PYBIND11_MODULE(libdescriptor, m) {
    m.doc() = "Python interface to libdescriptor";

    py::enum_<AvailableDescriptor>(m, "AvailableDescriptors")
            .value("Bispectrum", AvailableDescriptor::KindBispectrum)
            .value("SymmetryFunctions", AvailableDescriptor::KindSymmetryFunctions);

    py::class_<DescriptorKind, PyDescriptorKind>(m, "DescriptorKind")
            .def(py::init<>())
            .def("init_descriptor", py::overload_cast<AvailableDescriptor>(&DescriptorKind::initDescriptor))
            .def("init_descriptor",
                 py::overload_cast<std::string &, AvailableDescriptor>(&DescriptorKind::initDescriptor))
            .def("compute", &DescriptorKind::compute)
            .def_readwrite("param_file", &DescriptorKind::descriptor_param_file)
            .def_readwrite("kind", &DescriptorKind::descriptor_kind)
            .def_readwrite("width", &DescriptorKind::width);

    m.def("compute", &compute, "Compute descriptor for complete configuration.");
    m.def("compute_single_atom", &compute_single_atom, "Compute descriptor for single atom of configuration.");
    m.def("gradient", &gradient, "Compute gradient of descriptor for complete configuration.");
    m.def("gradient_single_atom", &gradient_single_atom,
          "Compute gradient of descriptor for single atom configuration.");
}


