#ifndef BISPECTRUM_HPP_
#define BISPECTRUM_HPP_

#include "helper.hpp"
#include "Descriptors.hpp"
#include <vector>

using namespace Descriptor;

/*! \class BISPECTRUM_LOOPINDICES
 * \brief The structure for the Bispectrum loop indices
 *
 */
struct BISPECTRUM_LOOPINDICES {
    int j1;
    int j2;
    int j;
};

/*!
 * \brief
 * This implementation is based on the method outlined
 * in Bartok[1], using formulae from VMK[2].
 *
 * For the Clebsch-Gordan coefficients, we convert the VMK half-integral
 * labels \f$ a, b, c, \alpha, \beta, \gamma \f$ to array offsets \f$ j_1, j_2, j, m_1, m_2, m \f$
 * using the following relations:
 *
 * $$ j_1 = 2*a $$
 * $$ j_2 = 2*b $$
 * $$ j =  2*c $$
 * $$ m_1 = \alpha+a $$ $$  2*\alpha = 2*m_1 - j_1 $$
 * $$ m_2 = \beta+b  $$ $$  2*\beta = 2*m2 - j_2  $$
 * $$ m =  \gamma+c  $$ $$  2*\gamma = 2*m - j   $$
 *
 * in this way:
 *
 * $$ -a \le \alpha \le a $$
 * $$ -b \le \beta \le b  $$
 * $$ -c \le \gamma \le c $$
 *
 * becomes:
 *
 *$$ 0 \le m_1 \le j_1 $$
 *$$ 0 \le m_2 \le j_2 $$
 *$$ 0 \le m \le j   $$
 *
 * and the requirement that
 * a+b+c be integral implies that
 * j1+j2+j must be even.
 * The requirement that:
 *
 * gamma = alpha+beta
 *
 * becomes:
 *
 * 2*m - j = 2*m1 - j1 + 2*m2 - j2
 *
 * Similarly, for the Wigner U-functions U(J,m,m') we
 * convert the half-integral labels J,m,m' to
 * array offsets j,ma,mb:
 *
 * j = 2*J
 * ma = J+m
 * mb = J+m'
 *
 * so that:
 *
 * 0 <= j <= 2*Jmax
 * 0 <= ma, mb <= j.
 *
 * For the bispectrum components B(J1,J2,J) we convert to:
 *
 * j1 = 2*J1
 * j2 = 2*J2
 * j = 2*J
 *
 * and the requirement:
 *
 * |J1-J2| <= J <= J1+J2, for j1+j2+j integral
 *
 * becomes:
 *
 * |j1-j2| <= j <= j1+j2, for j1+j2+j even integer
 *
 * or
 *
 * j = |j1-j2|, |j1-j2|+2,...,j1+j2-2,j1+j2
 *
 * [1] Albert Bartok-Partay, "Gaussian Approximation..."
 * Doctoral Thesis, Cambrindge University, (2009)
 * [2] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii,
 * "Quantum Theory of Angular Momentum," World Scientific (1988)
 *
 */
class Bispectrum final: public DescriptorKind {
    // final because otherwise I would need "virtual" destructor on DescriptorKind, and that makes enzyme spit out
    // warnings
public:
    Bispectrum(std::string &file_name);

    void initFromFile(std::string &file_name);

    Bispectrum() = default;

    /*!
   * \brief Construct a new Bispectrum object
   *
   * \param rfac0_in
   * \param twojmax_in
   * \param diagonalstyle_in
   * \param use_shared_arrays_in
   * \param rmin0_in
   * \param switch_flag_in
   * \param bzero_flag_in
   */
    Bispectrum(double rfac0_in,
               int twojmax_in,
               int diagonalstyle_in,
               int use_shared_arrays_in,
               double rmin0_in,
               int switch_flag_in,
               int bzero_flag_in);


    void build_indexlist();

    void init();

    double memory_usage();

    /*!
      * \brief Computes bispectrum for a set of atoms
      *
      * For example eq(5) of ``Gaussian Approximation Potentials: The Accuracy of
      * Quantum Mechanics, without the Electrons``, by Gabor Csany
      *
      * \param index
      * \param n_atoms
      * \param species
      * \param neigh_list
      * \param number_of_neigh
      * \param coords
      * \param desc
      */
    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neigh_list,
                 int number_of_neigh,
                 double_vector &coords,
                 double_vector &desc) override;

    void set_cutoff(const char *name, std::size_t Nspecies, double const *rcuts_in);

    void set_weight(int Nspecies, double const *weight_in);

//  void set_radius(int const Nspecies, double const * radius_in);
    void compute_ui(int jnum);

    void compute_zi();

    void compute_bi();

    void copy_bi2bvec();

//  void compute_duidrj(double const * rij_in,
//                      double const wj_in,
//                      double const rcut_in);
//  void compute_dbidrj();
//  void compute_dbidrj_nonsymm();
//  void copy_dbi2dbvec();
    double_scalar compute_sfac(double_scalar r, double_scalar rcut_in);

    double_scalar compute_dsfac(double_scalar r, double_scalar rcut_in);

    void grow_rij(int newnmax);

    void clone_empty(DescriptorKind *descriptorKind);

    void set_species(int n_species_) { n_species = n_species_; }

    int get_width();

private:
    inline double_scalar factorial(int n);

    void create_twojmax_arrays();

    void init_clebsch_gordan();

    void init_rootpqarray();

    inline void jtostr(char *str_out, int j);

    inline void mtostr(char *str_out, int j, int m);

    void print_clebsch_gordan(FILE *file);

    void zero_uarraytot();

    void addself_uarraytot(double_scalar wself_in);

    void add_uarraytot(double_scalar r, double_scalar wj_in, double_scalar rcut_in);

    void compute_uarray(double_scalar x, double_scalar y, double_scalar z, double_scalar z0, double_scalar r);

    inline double_scalar deltacg(int j1, int j2, int j);

    int compute_ncoeff();

public:
    int ncoeff;
    std::vector<double_scalar> bvec;
    Array2D<double_scalar> dbvec;
    Array2D<double_scalar> rij;
    std::vector<int> inside;
    std::vector<double_scalar> wj;
    std::vector<double_scalar> rcutij;
    int nmax;
    int twojmax;
    int diagonalstyle;
    Array3D<double_scalar> uarraytot_r;
    Array3D<double_scalar> uarraytot_i;
    Array5D<double_scalar> zarray_r;
    Array5D<double_scalar> zarray_i;

//  Array3D<double_scalar> uarraytot_r_b;
//  Array3D<double_scalar> uarraytot_i_b;

//  Array5D<double_scalar> zarray_r_b;
//  Array5D<double_scalar> zarray_i_b;

    Array3D<double_scalar> uarray_r;
    Array3D<double_scalar> uarray_i;

private:
    int n_species;
    Array2D<double_scalar> rcuts;
    std::vector<double_scalar> wjelem;
    double_scalar rmin0;
    double_scalar rfac0;
    std::vector<BISPECTRUM_LOOPINDICES> idxj;
    int idxj_max;
    Array5D<double_scalar> cgarray;
    std::vector<double_scalar> rootpqarray;
    Array3D<double_scalar> barray;
    Array4D<double_scalar> duarray_r;
    Array4D<double_scalar> duarray_i;
    Array4D<double_scalar> dbarray;
    static const int nmaxfactorial = 167;
    static double_scalar const nfac_table[];
    int use_shared_arrays;
    int switch_flag;
    double_scalar wself;
    int bzero_flag;
    std::vector<double_scalar> bzero;
};

#endif  // BISPECTRUM_HPP_
