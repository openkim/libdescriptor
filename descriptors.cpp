#include "SymFun.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

namespace py = pybind11;

PYBIND11_MODULE(descriptors, m) {
    m.doc() = "Symmetry function descriptor for ANN potential.";

    py::class_<SymmetryFunctionParams>(m, "SymmetryFunctionParams")
            .def(py::init<>())
//      .def(py::init<SymmetryFunctionParams const &>())

            .def("get_num_descriptors", &SymmetryFunctionParams::get_num_descriptors)

            .def(
                    "set_cutoff",
                    [](SymmetryFunctionParams &d, char *name, py::array_t<double> rcuts) {
                        d.set_cutoff(name, rcuts.shape(0), rcuts.data(0));
                        return;
                    },
                    py::arg("name"),
                    py::arg("rcuts").noconvert())

            .def(
                    "add_descriptor",
                    [](SymmetryFunctionParams &d, char *name, py::array_t<double> values) {
                        auto rows = values.shape(0);
                        auto cols = values.shape(1);
                        d.add_descriptor(name, values.data(0), rows, cols);
                        return;
                    },
                    py::arg("name"),
                    py::arg("values").noconvert())
            .def(
                    "num_params",
                    [](SymmetryFunctionParams &d) {
                        return py::array(py::buffer_info(
                                d.num_params_.data(),
                                sizeof(int),
                                py::format_descriptor<int>::format(),
                                1,
                                {d.num_params_.size()},
                                {sizeof(int)}
                        ));
                    }
            );

    m.def(
            "symmetry_function_atomic",
            [](
                    int i,
                    py::array_t<double> coords,
                    py::array_t<int> particleSpecies,
                    py::array_t<int> neighlist,
                    int width,
                    SymmetryFunctionParams &symparm
            ) {
                int numnei = neighlist.shape(0);
                std::vector<double> zeta(width, 0.0);
                symmetry_function_atomic(
                        i,
                        coords.data(0),
                        particleSpecies.data(0),
                        neighlist.data(0),
                        numnei,
                        zeta.data(),
                        &symparm);
                auto zeta_py = py::array(py::buffer_info(
                        zeta.data(),
                        sizeof(double),
                        py::format_descriptor<double>::format(),
                        1,
                        {width},
                        {sizeof(double)}
                ));
                return zeta_py;
            },
            py::arg("i"),
            py::arg("coords").noconvert(),
            py::arg("particleSpecies").noconvert(),
            py::arg("neighlist").noconvert(),
            py::arg("width"),
            py::arg("symparm"),
            "Return zeta");
    m.def(
            "grad_symmetry_function_atomic",
            [](
                    int i,
                    py::array_t<double> coords,
                    py::array_t<int> particleSpecies,
                    py::array_t<int> neighlist,
                    int width,
                    SymmetryFunctionParams &symparm,
                    py::array_t<double> dE_dzeta
            ) {
                int numnei = neighlist.shape(0);
                std::vector<double> d_coords(coords.size(), 0.0);
                std::vector<double> zeta(width, 0.0);
                grad_symmetry_function_atomic(i,
                                              coords.data(0),
                                              d_coords.data(),
                                              particleSpecies.data(0),
                                              neighlist.data(0),
                                              numnei,
                                              zeta.data(),
                                              dE_dzeta.data(0),
                                              &symparm);
                auto dr_dzeta = py::array(py::buffer_info(
                        d_coords.data(),
                        sizeof(double),
                        py::format_descriptor<double>::format(),
                        1,
                        {static_cast<int>(coords.size())},
                        {sizeof(double)}
                ));
                return dr_dzeta;
            },
            py::arg("i"),
            py::arg("coords").noconvert(),
            py::arg("particleSpecies").noconvert(),
            py::arg("neighlist").noconvert(),
            py::arg("width"),
            py::arg("symparm"),
            py::arg("dE_dzeta").noconvert(),
            "Return derivative"
    );
}
