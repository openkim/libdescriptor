#ifndef KLIFF_DESCRIPTOR_HPP_
#define KLIFF_DESCRIPTOR_HPP_

#include "helper.hpp"

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


class SymmetryFunctionParams {
public:
    SymmetryFunctionParams();

    ~SymmetryFunctionParams();

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
    bool has_three_body_;
};

#undef DIM

inline double cut_cos(double r, double rcut);


void symmetry_function_atomic(int i,
                              double const *coords,
                              int const *particleSpeciesCodes,
                              int const *neighlist,
                              int numnei,
                              double * desc,
                              SymmetryFunctionParams *SymParam);

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

void __enzyme_autodiff(void (*)(
        int const,
        double const *,
        int const *,
        int const *,
        int const,
        double *const,
        SymmetryFunctionParams *),
                       int, int,
                       int, double const *, double const *,
                       int, int const *,
                       int, int const *,
                       int, int,
                       int, double *, double const *,
                       int, SymmetryFunctionParams *);

void grad_symmetry_function_atomic(int  i,
                                   double const *coords,
                                   double const *d_coords,
                                   int const *particleSpeciesCodes,
                                   int const *neighlist,
                                   int  numnei,
                                   double * desc,
                                   double const *d_grad_loss_zeta,
                                   SymmetryFunctionParams *SymParam);

#endif  // KLIFF_DESCRIPTOR_HPP_
