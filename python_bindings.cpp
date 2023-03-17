#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <vector>
#include <string>
#include <iostream>
#include "Descriptors.hpp"

namespace py = pybind11;
using namespace Descriptor;

// Trampoline class for DescriptorKind pure virtual compute function
class PyDescriptorKind : public DescriptorKind {
public:
    using DescriptorKind::DescriptorKind;

    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neighbor_lists,
                 int number_of_neighbors,
                 double_vector &coordinates,
                 double_vector &desc) override {
        PYBIND11_OVERRIDE_PURE(void, DescriptorKind, compute,
                               index, n_atoms, species, neighbor_lists, number_of_neighbors, coordinates, desc
        );
    };
};


PYBIND11_MODULE(libdescriptor, m) {
    m.doc() = "Python interface to libdescriptor";

    // Python bindings for AvailableDescriptor kinds.
    py::enum_<AvailableDescriptor>(m, "AvailableDescriptors")
            .value("Bispectrum", AvailableDescriptor::KindBispectrum)
            .value("SymmetryFunctions", AvailableDescriptor::KindSymmetryFunctions);

    // Python bindings for DescriptorKind
    py::class_<DescriptorKind, PyDescriptorKind>(m, "DescriptorKind")
            .def(py::init<>()) // Default constructor

                    // Return default constructor, mostly needed for creating empty objects, may need in future for
                    // hyperparameter optimization
            .def("init_descriptor", py::overload_cast<AvailableDescriptor>(&DescriptorKind::initDescriptor))

                    // Return descriptor object created from reading the input file
            .def("init_descriptor",
                 py::overload_cast<std::string &, AvailableDescriptor>(&DescriptorKind::initDescriptor))

                    // Return descriptor object created from Symmetry Functions, using direct parameters
            .def("init_descriptor", py::overload_cast<AvailableDescriptor,
                    std::vector<std::string> *, std::string *, double *, std::vector<std::string> *, std::vector<int> *,
                    std::vector<double> *>(&DescriptorKind::initDescriptor))

//            .def("dummy", [](AvailableDescriptor a,std::vector<std::string> *b, std::string *c, double *d, std::vector<std::string> *e, std::vector<int> *f,
//                    std::vector<double> *g){return DescriptorKind::initDescriptor(a,b,c,d,e,f,g);})

                    // Return descriptor object created from Bispectrum, using direct parameters
            .def("init_descriptor", py::overload_cast<AvailableDescriptor, double, int, int, int, double, int,
                 int, double *, std::vector<std::string> *, std::vector<double> *>(&DescriptorKind::initDescriptor))

                    // Compute function for calculating the descriptor
            .def("compute",
                 [](DescriptorKind &ds, int index, py::array_t<int, py::array::c_style | py::array::forcecast> &species,
                    py::array_t<int, py::array::c_style | py::array::forcecast> &neighbors,
                    py::array_t<double, py::array::c_style | py::array::forcecast> &coords) {
                     int n_atoms = static_cast<int>(coords.shape(0));
                     int n_neigh = static_cast<int>(neighbors.shape(0));
                     double_vector desc(ds.width);
                     for (int i = 0; i < ds.width; i++) desc[i] = 0.0;
                     double_vector coords_vec(n_atoms * 3);
                        for (int i = 0; i < n_atoms; i++) {
                            coords_vec[i * 3] = coords.at(i, 0);
                            coords_vec[i * 3 + 1] = coords.at(i, 1);
                            coords_vec[i * 3 + 2] = coords.at(i, 2);
                        }
                     ds.compute(index,
                                n_atoms,
                                const_cast<int *>(species.data(0)),
                                const_cast<int *>(neighbors.data(0)),
                                n_neigh,
                                coords_vec,
                                desc);

                     // copy descriptor to numpy array
                     py::array_t<double> desc_array(ds.width);
                     for (int i = 0; i < ds.width; i++) {
                            desc_array.mutable_at(i) = static_cast<double>(desc[i]);
                     }
                     return desc_array;
                 }, py::return_value_policy::take_ownership)
            .def_readwrite("param_file", &DescriptorKind::descriptor_param_file)
            .def_readwrite("kind", &DescriptorKind::descriptor_kind)
            .def_readwrite("width", &DescriptorKind::width);


    // m.def("compute", &compute, "Compute descriptor for complete configuration.");
    // Calculate descriptor for a single atom using the descriptor supplied in DescriptorKind argument
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

    // Calculate gradient of descriptor for a single atom using the descriptor supplied in DescriptorKind argument
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

    // Calculate gradient for a single atom using the numerical differentiation
    m.def("num_gradient_single_atom",
          [](DescriptorKind &ds, int index, int n_contributing_atoms, py::array_t<int, py::array::c_style | py::array::forcecast> &species,
             py::array_t<int, py::array::c_style | py::array::forcecast> &neighbors,
             py::array_t<double, py::array::c_style | py::array::forcecast> &coordinates,
             py::array_t<double, py::array::c_style | py::array::forcecast> &dE_ddesc) {
              int n_total_atoms = static_cast<int>(coordinates.shape(0));
              int n_neigh = static_cast<int>(neighbors.shape(0));
              auto d_coordinates = new double[3];
              for (int i = 0; i < 3; i++) d_coordinates[i] = 0.0;
              num_gradient_single_atom(index,
                                       n_contributing_atoms,
                                       n_total_atoms,
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
