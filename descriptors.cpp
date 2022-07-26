#include "descriptors.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <fstream>
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


Descriptor::Descriptor() {
    descriptor_map.insert(std::make_pair("SymFun",reinterpret_cast<void *>(new SymmetryFunctionParams)));
}

Descriptor::Descriptor(std::string& descriptor_name) {
    descriptor_kind = descriptor_name;
    descriptor_map.insert(std::make_pair("SymFun",reinterpret_cast<void *>(new SymmetryFunctionParams)));
}

Descriptor::Descriptor(std::string& descriptor_name, std::string& descriptor_params) {
    descriptor_kind = descriptor_name;
    descriptor_param_file = descriptor_name;
    descriptor_map.insert(std::make_pair("SymFun",reinterpret_cast<void *>(new SymmetryFunctionParams)));
    initDescriptor();
}

void Descriptor::initDescriptor() {
    if (descriptor_kind == "SymFun"){
        std::fstream file_ptr(descriptor_param_file);
        std::string placeholder_string;
        int n_species;

        // Ignore comments
        do {
            std::getline(file_ptr, placeholder_string);
        } while (placeholder_string[0] == '#');
        n_species = std::stoi(placeholder_string);

        // blank line
        std::getline(file_ptr, placeholder_string);
        // Ignore comments
        do {
            std::getline(file_ptr, placeholder_string);
        } while (placeholder_string[0] == '#');

        double * cutoff_matrix = new double[n_species * n_species];
        for (int i=0;i<n_species;i++){
            std::getline(file_ptr, placeholder_string);
            for (int j=0;j<n_species;j++){
                auto pos = placeholder_string.find(' ');
                cutoff_matrix[n_species * i + j] = stod(placeholder_string.substr(0, pos));
                if (pos != std::string::npos) placeholder_string.erase(0, pos + 1);
            }
        }

        // blank line
        std::getline(file_ptr, placeholder_string);
        // Ignore comments
        do {
            std::getline(file_ptr, placeholder_string);
        } while (placeholder_string[0] == '#');
        std::string cutoff_function = placeholder_string;

        // blank line
        std::getline(file_ptr, placeholder_string);
        // Ignore comments
        do {
            std::getline(file_ptr, placeholder_string);
        } while (placeholder_string[0] == '#');
        int width = std::stoi(placeholder_string);

        // blank line
        std::getline(file_ptr, placeholder_string);
        // Ignore comments
        do {
            std::getline(file_ptr, placeholder_string);
        } while (placeholder_string[0] == '#');
        int n_func = std::stoi(placeholder_string);

        std::vector<std::string> sym_func_list;
        for (int i = 0; i < n_func; i++) {
            std::getline(file_ptr, placeholder_string);
            sym_func_list.push_back(placeholder_string);
        }

        std::vector<int> sym_func_lengths;
        for (int i = 0; i < n_func; i++) {
            std::getline(file_ptr, placeholder_string);
            sym_func_lengths.push_back(std::stoi(placeholder_string));
        }

        std::vector<std::vector<double>> sym_func_elements;
        for (int i = 0; i < n_func; i++) {
            //blank line
            std::getline(file_ptr, placeholder_string);
            std::vector<double> tmp_desc_list;
            for (int j = 0; j < sym_func_lengths[i]; j++){
                std::getline(file_ptr, placeholder_string);
                tmp_desc_list.push_back(std::stod(placeholder_string));
            }
            sym_func_elements.push_back(std::move(tmp_desc_list));
        }

        for (int i = 0; i < n_func ; i++){
            if (sym_func_list[i] == "g2"){
                for (int j = 0; j < sym_func_elements[i].size(); j = j+2) sym_func_elements[i][j] /= (bhor2ang * bhor2ang);
            } else if (sym_func_list[i] == "g4") {
                for (int j = 2; j < sym_func_elements[i].size(); j = j+3) sym_func_elements[i][j] /= (bhor2ang * bhor2ang);
            }
        }

        // blank line
        std::getline(file_ptr, placeholder_string);
        // Ignore comments
        do {
            std::getline(file_ptr, placeholder_string);
        } while (placeholder_string[0] == '#');

        std::vector<int> dims;
        for (int i=0; i < n_func * 2; i++){
            dims.push_back(std::stoi(placeholder_string));
            std::getline(file_ptr, placeholder_string);
        }

        auto sf = reinterpret_cast<SymmetryFunctionParams *> (descriptor_map["SymFun"]);
        sf->set_cutoff(cutoff_function.c_str(), 1, cutoff_matrix);
        for (int i =0 ; i< n_func; i++){
            sf->add_descriptor(sym_func_list[i].c_str(), sym_func_elements[i].data(), dims[2 *i], dims[2*i+1] );
        }
        sf->width = width;
        delete[] cutoff_matrix;
    }
}

