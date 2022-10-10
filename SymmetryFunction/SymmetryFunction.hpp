#ifndef SYMMETRY_FUNCTION_HPP_
#define SYMMETRY_FUNCTION_HPP_

#include "helper.hpp"
#include "descriptors.hpp"

#include <numeric>

#include <cmath>
#include <cstring>
#include <map>
#include <string>

#ifdef DIM
#undef DIM
#endif

#define DIM 3

typedef double VectorOfSizeDIM[DIM];

using namespace Descriptor;

class SymmetryFunctions : public DescriptorKind {
public:
    SymmetryFunctions(std::string &);
    SymmetryFunctions()=default;

    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neigh_list,
                 int number_of_neigh,
                 double *coords,
                 double *zeta) override;

    void set_length() override {length = 51;};

private:

    bool has_three_body_;

    inline void set_species(std::vector<std::string> &species);

    inline void get_species(std::vector<std::string> &species);

    inline int get_num_species();

    int width;

    void set_cutoff(char const *name,
                    std::size_t Nspecies,
                    double const *rcut_2D);

    inline double get_cutoff(int iCode, int jCode);

    void add_descriptor(char const *name,
                        double const *values,
                        int row,
                        int col);

    int get_num_descriptors();

    std::vector<std::string> species_;
    std::vector<int> name_;
    std::vector<int> starting_index_;
    Array2D<double> rcut_2D_;
    std::vector<Array2D<double> > params_;
    std::vector<int> num_param_sets_;
    std::vector<int> num_params_;
};

inline double cut_cos(double r, double rcut);

void sym_g1(double r, double rcut, double &phi);

void sym_g2(double eta,
            double Rs,
            double r,
            double rcut,
            double &phi);

void sym_g3(double kappa, double r, double rcut, double &phi);


void sym_g4(double zeta,
            double lambda,
            double eta,
            double const *r,
            double const *rcut,
            double &phi);

void sym_g5(double zeta,
            double lambda,
            double eta,
            double const *r,
            double const *rcut,
            double &phi);

#endif  // SYMMETRY_FUNCTION_HPP_
