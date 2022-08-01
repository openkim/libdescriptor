#include "descriptors.hpp"
#include <iostream>
#include <fstream>
#include <vector>

Descriptor::Descriptor() {
    descriptor_map.insert(std::make_pair("SymmetryFunction",reinterpret_cast<void *>(new SymmetryFunctionParams)));
}

Descriptor::Descriptor(std::string& descriptor_name) {
    descriptor_kind = descriptor_name;
    descriptor_map.insert(std::make_pair("SymmetryFunction",reinterpret_cast<void *>(new SymmetryFunctionParams)));
}

Descriptor::Descriptor(std::string& descriptor_name, std::string& descriptor_params) {
    descriptor_kind = descriptor_name;
    descriptor_param_file = descriptor_params;
    descriptor_map.insert(std::make_pair("SymmetryFunction",reinterpret_cast<void *>(new SymmetryFunctionParams)));
    initDescriptor();
}

void Descriptor::initDescriptor() {

    if (descriptor_kind == "SymmetryFunction"){
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
            for (int j=0;j<n_species;j++){
                auto pos = placeholder_string.find(' ');
                *( cutoff_matrix + n_species * i + j) = std::stod(placeholder_string.substr(0, pos));
                if (pos != std::string::npos) placeholder_string.erase(0, pos + 1);
            }
            std::getline(file_ptr, placeholder_string);
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

        auto sf = reinterpret_cast<SymmetryFunctionParams *> (descriptor_map["SymmetryFunction"]);

        sf->set_cutoff(cutoff_function.c_str(), n_species, cutoff_matrix);

        for (int i =0 ; i< n_func; i++){
            sf->add_descriptor(sym_func_list[i].c_str(), sym_func_elements[i].data(), dims[2 *i], dims[2*i+1] );
        }
        sf->width = width;
        delete[] cutoff_matrix;
    }
}

Descriptor::~Descriptor() {};