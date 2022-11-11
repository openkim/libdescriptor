#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>
#include <string>
#include <iostream>
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
            .def("compute",
                 [](DescriptorKind &ds, int index, py::array_t<int, py::array::c_style | py::array::forcecast> &species,
                    py::array_t<int, py::array::c_style | py::array::forcecast> &neighbors,
                    py::array_t<double, py::array::c_style | py::array::forcecast> &coords) {
                     int n_atoms = static_cast<int>(coords.shape(0));
                     int n_neigh = static_cast<int>(neighbors.shape(0));
                     auto desc = new double[ds.width];
                     for (int i = 0; i < ds.width; i++) desc[i] = 0.0;
                     ds.compute(index,
                                n_atoms,
                                const_cast<int *>(species.data(0)),
                                const_cast<int *>(neighbors.data(0)),
                                n_neigh,
                                const_cast<double *>(coords.data(0)),
                                desc);

                     py::array_t<double> desc_array(ds.width, desc);
                     return desc_array;
                 }, py::return_value_policy::take_ownership)
            .def_readwrite("param_file", &DescriptorKind::descriptor_param_file)
            .def_readwrite("kind", &DescriptorKind::descriptor_kind)
            .def_readwrite("width", &DescriptorKind::width);

    //#TODO add individual descriptors

    // m.def("compute", &compute, "Compute descriptor for complete configuration.");
    m.def("compute_single_atom",
          [](DescriptorKind &ds, int index, py::array_t<int, py::array::c_style | py::array::forcecast> &species,
             py::array_t<int, py::array::c_style | py::array::forcecast> &neighbors,
             py::array_t<double, py::array::c_style | py::array::forcecast> &coordinates) {
              int n_atoms = static_cast<int>(coordinates.shape(0));
              int n_neigh = static_cast<int>(neighbors.shape(0));
              auto desc = new double[ds.width];
              for (int i = 0; i < ds.width; i++) desc[i] = 0.0;
              compute_single_atom(index,
                                  n_atoms,
                                  const_cast<int *>(species.data(0)),
                                  const_cast<int *>(neighbors.data(0)),
                                  n_neigh,
                                  const_cast<double *>(coordinates.data(0)),
                                  desc,
                                  &ds);

              py::array_t<double> desc_array(ds.width, desc);
              return desc_array;
          }, py::return_value_policy::take_ownership,
          "Compute descriptor for single atom of configuration.");
    // m.def("gradient", &gradient, "Compute gradient of descriptor for complete configuration.");
    m.def("gradient_single_atom",
          [](DescriptorKind &ds, int index, py::array_t<int, py::array::c_style | py::array::forcecast> &species,
             py::array_t<int, py::array::c_style | py::array::forcecast> &neighbors,
             py::array_t<double, py::array::c_style | py::array::forcecast> &coordinates,
             py::array_t<double, py::array::c_style | py::array::forcecast> &desc,
             py::array_t<double, py::array::c_style | py::array::forcecast> &dE_ddesc) {
              int n_atoms = static_cast<int>(coordinates.shape(0));
              int n_neigh = static_cast<int>(neighbors.shape(0));
              auto d_coordinates = new double[coordinates.size()];
              for (int i = 0; i < coordinates.size(); i++) d_coordinates[i] = 0.0;
              gradient_single_atom(index,
                                   n_atoms,
                                   const_cast<int *>(species.data(0)),
                                   const_cast<int *>(neighbors.data(0)),
                                   n_neigh,
                                   const_cast<double *>(coordinates.data(0)),
                                   d_coordinates,
                                   const_cast<double *>(desc.data(0)),
                                   const_cast<double *>(dE_ddesc.data(0)),
                                   &ds);
              py::array_t<double> d_coord_array({coordinates.shape(0), coordinates.shape(1)}, d_coordinates);
              return d_coord_array;
          }, py::return_value_policy::take_ownership,
          "Compute gradient of descriptor for single atom configuration.");
    // TODO generate_one_atom, which is compatible with current kliff
    m.def("num_gradient_single_atom",
          [](DescriptorKind &ds, int index, py::array_t<int, py::array::c_style | py::array::forcecast> &species,
             py::array_t<int, py::array::c_style | py::array::forcecast> &neighbors,
             py::array_t<double, py::array::c_style | py::array::forcecast> &coordinates,
             py::array_t<double, py::array::c_style | py::array::forcecast> &dE_ddesc) {
              int n_atoms = static_cast<int>(coordinates.shape(0));
              int n_neigh = static_cast<int>(neighbors.shape(0));
              auto d_coordinates = new double[3];
              for (int i = 0; i < 3; i++) d_coordinates[i] = 0.0;
              num_gradient_single_atom(index,
                                   n_atoms,
                                   const_cast<int *>(species.data(0)),
                                   const_cast<int *>(neighbors.data(0)),
                                   n_neigh,
                                   const_cast<double *>(coordinates.data(0)),
                                   d_coordinates,
                                   const_cast<double *>(dE_ddesc.data(0)),
                                   &ds);
              py::array_t<double> d_coord_array(3, d_coordinates);
              return d_coord_array;
          }, py::return_value_policy::take_ownership,
          "Compute gradient of descriptor for single atom configuration.");
}
