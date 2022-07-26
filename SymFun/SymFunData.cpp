#include "SymFun.hpp"
#include <cstring>
#include <iostream>

#ifdef MY_PI
#undef MY_PI
#endif

#define MY_PI 3.1415926535897932

#ifdef DIM
#undef DIM
#endif

#define DIM 3

#ifdef MAXLINE
#undef MAXLINE
#endif

#define MAXLINE 20480

#ifdef LOG_ERROR
#undef LOG_ERROR
#endif

#define LOG_ERROR(msg)                                           \
  {                                                              \
    std::ostringstream ss;                                       \
    ss << msg;                                                   \
    std::string _Messagef_(FormatMessageFileLineFunctionMessage( \
        "Error ", __FILE__, __LINE__, __FUNCTION__, ss.str()));  \
    std::cerr << _Messagef_;                                     \
  }


SymmetryFunctionParams::SymmetryFunctionParams() : has_three_body_(false) {
}

SymmetryFunctionParams::~SymmetryFunctionParams() {}

inline void SymmetryFunctionParams::set_species(std::vector<std::string> &species) {
    species_.resize(species.size());
    std::copy(species.begin(), species.end(), species_.begin());
};

inline void SymmetryFunctionParams::get_species(std::vector<std::string> &species) {
    species.resize(species_.size());
    std::copy(species_.begin(), species_.end(), species.begin());
};

inline int SymmetryFunctionParams::get_num_species() { return species_.size(); }

void SymmetryFunctionParams::set_cutoff(char const *name,
                                        std::size_t const Nspecies,
                                        double const *rcut_2D) {
    (void) name;   // to avoid unused warning
    rcut_2D_.resize(Nspecies, Nspecies, rcut_2D);
}

inline double SymmetryFunctionParams::get_cutoff(int const iCode, int const jCode) {
    return rcut_2D_(iCode, jCode);
};

void SymmetryFunctionParams::add_descriptor(char const *name,
                                            double const *values,
                                            int const row,
                                            int const col) {
    if (strcmp(name, "g1") == 0) { name_.push_back(1); };
    if (strcmp(name, "g2") == 0) { name_.push_back(2); };
    if (strcmp(name, "g3") == 0) { name_.push_back(3); };
    if (strcmp(name, "g4") == 0) { name_.push_back(4); };
    if (strcmp(name, "g5") == 0) { name_.push_back(5); };

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

int SymmetryFunctionParams::get_num_descriptors() {
    return std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
}

#undef LOG_ERROR
#undef MAXLINE
#undef DIM
#undef MY_PI
