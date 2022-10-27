#ifndef SYMMETRY_FUNCTION_HPP_
#define SYMMETRY_FUNCTION_HPP_

#include "helper.hpp"
#include "Descriptors.hpp"

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
    // In KLIFF a utility will create in memory file for initiaization?
    // TODO create apropriate constructor
    explicit SymmetryFunctions(std::string &file_name);

    void initFromFile(std::string &file_name);

    SymmetryFunctions() {};

    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neigh_list,
                 int number_of_neigh,
                 double *coords,
                 double *zeta) override;

    void clone_empty(DescriptorKind *descriptorKind);

private:
    int n_species;

    bool has_three_body_;

    inline void set_species(std::vector<std::string> &species);

    inline void get_species(std::vector<std::string> &species);

    inline int get_num_species();

    void set_cutoff(char const *name,
                    std::size_t Nspecies,
                    double const *rcut_2D);

    inline double get_cutoff(int iCode, int jCode);

    void add_descriptor(char const *name,
                        double const *values,
                        int row,
                        int col);

    int get_num_descriptors();

    double bhor2ang = 0.529177;
    std::vector<std::string> species_;
    std::vector<int> name_;
    std::vector<int> starting_index_;
    Array2D<double> rcut_2D_;
    std::vector<Array2D<double> > params_;
    std::vector<int> num_param_sets_;
    std::vector<int> num_params_;
};

inline double cut_cos(double r, double rcut);

// Symmetry functions: Jorg Behler, J. Chem. Phys. 134, 074106, 2011.

/*!
 * \brief Radial symmetry function \c g1 suitable for describing the
 * radial environment of atom \c i
 *
 * \c g1 is the sum of the cutoff functions with respect to
 * all neighboring atoms
 *
 * \param r The distance between atoms \c i and \c j
 * \param rcut The cutoff radius
 * \param phi Radial symmetry function \c g1 value
 */
void sym_g1(double r, double rcut, double &phi);

/*!
 * \brief Radial symmetry function \c g2 suitable for describing the
 * radial environment of atom \c i
 *
 * \c g2 is the sum of Gaussians multiplied by cutoff functions. The shifted
 * \g2 is suitable to describe a spherical shell around the reference atom.
 *
 * \param eta Defines the width of the Gaussians
 * \param Rs Shifts the center of the Gaussians to a certain radial distance
 * \param r The distance between atoms \c i and \c j
 * \param rcut The cutoff radius
 * \param phi Radial symmetry function \c g2 value
 */
void sym_g2(double eta, double Rs, double r, double rcut, double &phi);

/*!
 * \brief Radial symmetry function \c g3 suitable for describing the
 * radial environment of atom \c i
 *
 * \c g3 is a damped cosine function with a period length adjusted
 * by parameter \c kappa
 *
 * \param kappa Adjusting the period length of the damped cosine function
 * \param r The distance between atoms \c i and \c j
 * \param rcut The cutoff radius
 * \param phi Radial symmetry function \c g3 value
 *
 * \note
 * Care must be taken when using the \c g3 function. Neighboring atoms can
 * cancel each other’s due to the existence of positive and negative function
 * values. It is recommended to use \c g3 in combination with other symmetry
 * functions.
 */
void sym_g3(double kappa, double r, double rcut, double &phi);

/*!
 * \brief Summations of cosine functions of the angles centered at atom \c i.
 *
 * The angular part must be symmetric with respect to \c 180 angle
 *
 * \param zeta Provides the angular resolution. High values yield a narrower
 * range of nonzero symmetry function values \param lambda Shifting the maxima
 * of the cosine function (can have the values +1 or -1) \param eta
 * Controlling the radial resolution of the radial function. \param r The
 * distance between atoms \c i and \c j \param rcut The cutoff radius \param
 * phi Function \c g4 value
 */
void sym_g4(double zeta, double lambda, double eta, double const *r, double const *rcut, double &phi);

/*!
 * \brief Summations of cosine functions of the angles centered at atom \c i.
 *
 * The angular part must be symmetric with respect to \c 180 angle.
 * In function \c g5 there is no constraint on the distance between atoms
 * resulting in a larger number of terms in the summation.
 *
 * \param zeta Provides the angular resolution. High values yield a narrower
 * range of nonzero symmetry function values \param lambda Shifting the maxima
 * of the cosine function (can have the values +1 or -1) \param eta
 * Controlling the radial resolution of the radial function. \param r The
 * distance between atoms \c i and \c j \param rcut The cutoff radius \param
 * phi Function \c g4 value
 */
void sym_g5(double zeta, double lambda, double eta, double const *r, double const *rcut, double &phi);

#endif  // SYMMETRY_FUNCTION_HPP_
