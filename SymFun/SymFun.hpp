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


//using CutoffFunction = double (*)(double, double);
class SymmetryFunctionParams {
public:
    SymmetryFunctionParams();

//  SymmetryFunctionParams(const SymmetryFunctionParams &sym_fun_param1);
//  SymmetryFunctionParams(int zero,const SymmetryFunctionParams &sym_fun_param1);
    ~SymmetryFunctionParams();

    inline void set_species(std::vector<std::string> &species);

    inline void get_species(std::vector<std::string> &species);

    inline int get_num_species();

    void set_cutoff(char const *name,
                    std::size_t const Nspecies,
                    double const *rcut_2D);

    inline double get_cutoff(int const iCode, int const jCode);

    void add_descriptor(char const *name,
                        double const *values,
                        int const row,
                        int const col);

    int get_num_descriptors();

    void echo_input();

//  CutoffFunction cutoff_func_;

    std::vector<std::string> species_;
    std::vector<int> name_;
    std::vector<int> starting_index_;
    Array2D<double> rcut_2D_;
    std::vector<Array2D<double> > params_;
    std::vector<int> num_param_sets_;
    std::vector<int> num_params_;
//  std::vector<double> feature_mean_;
//  std::vector<double> feature_std_;
    bool has_three_body_;
//  bool normalize_;
};

#undef DIM

inline double cut_cos(double const r, double const rcut);


void symmetry_function_atomic(int const i,
                              double const *coords,
                              int const *particleSpeciesCodes,
                              int const *neighlist,
                              int const numnei,
                              double *const desc,
                              SymmetryFunctionParams *SymParam);

void sym_g1(double const r, double const rcut, double &phi);

void sym_g2(double const eta,
            double const Rs,
            double const r,
            double const rcut,
            double &phi);

void sym_g3(double const kappa, double const r, double const rcut, double &phi);


void sym_g4(double const zeta,
            double const lambda,
            double const eta,
            double const *r,
            double const *rcut,
            double &phi);

void sym_g5(double const zeta,
            double const lambda,
            double const eta,
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
                       int, int const,
                       int, double const *, double const *,
                       int, int const *,
                       int, int const *,
                       int, int const,
                       int, double *const, double const *,
                       int, SymmetryFunctionParams *);

void grad_symmetry_function_atomic(int const i,
                                   double const *coords,
                                   double const *d_coords,
                                   int const *particleSpeciesCodes,
                                   int const *neighlist,
                                   int const numnei,
                                   double *const desc,
                                   double const *d_grad_loss_zeta,
                                   SymmetryFunctionParams *SymParam);

void test();

#endif  // KLIFF_DESCRIPTOR_HPP_
