#ifndef BISPECTRUM_HPP_
#define BISPECTRUM_HPP_

#include "helper.hpp"
#include "descriptors.hpp"

using namespace Descriptor;

struct BISPECTRUM_LOOPINDICES {
    int j1;
    int j2;
    int j;
};

class Bispectrum : public DescriptorKind {
public:
    Bispectrum(std::string &);
    Bispectrum() = default;
    Bispectrum(double rfac0_in,
               int twojmax_in,
               int diagonalstyle_in,
               int use_shared_arrays_in,
               double rmin0_in,
               int switch_flag_in,
               int bzero_flag_in);

    ~Bispectrum();

    void build_indexlist();

    void init();

    double memory_usage();

    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neigh_list,
                 int number_of_neigh,
                 double *coords,
                 double *desc) override;

    void set_cutoff(char *name, std::size_t Nspecies, double const *rcuts_in);

    void set_length() override {length = 51;}

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
    double compute_sfac(double r, double rcut_in);

    double compute_dsfac(double r, double rcut_in);

    void grow_rij(int newnmax);

private:
    inline double factorial(int n);

    void create_twojmax_arrays();

    void init_clebsch_gordan();

    void init_rootpqarray();

    inline void jtostr(char *str_out, int j);

    inline void mtostr(char *str_out, int j, int m);

    void print_clebsch_gordan(FILE *file);

    void zero_uarraytot();

    void addself_uarraytot(double wself_in);

    void add_uarraytot(double r, double wj_in, double rcut_in);

    void compute_uarray(double x,
                        double y,
                        double z,
                        double z0,
                        double r);

    inline double deltacg(int j1, int j2, int j);

    int compute_ncoeff();

public:
    int ncoeff;
    std::vector<double> bvec;
    Array2D<double> dbvec;
    Array2D<double> rij;
    std::vector<int> inside;
    std::vector<double> wj;
    std::vector<double> rcutij;
    int nmax;
    int twojmax;
    int diagonalstyle;
    Array3D<double> uarraytot_r;
    Array3D<double> uarraytot_i;
    Array5D<double> zarray_r;
    Array5D<double> zarray_i;

//  Array3D<double> uarraytot_r_b;
//  Array3D<double> uarraytot_i_b;

//  Array5D<double> zarray_r_b;
//  Array5D<double> zarray_i_b;

    Array3D<double> uarray_r;
    Array3D<double> uarray_i;

private:
    Array2D<double> rcuts;
    std::vector<double> wjelem;
    double rmin0;
    double rfac0;
    std::vector<BISPECTRUM_LOOPINDICES> idxj;
    int idxj_max;
    Array5D<double> cgarray;
    Array2D<double> rootpqarray;
    Array3D<double> barray;
    Array4D<double> duarray_r;
    Array4D<double> duarray_i;
    Array4D<double> dbarray;
    static const int nmaxfactorial = 167;
    static double const nfac_table[];
    int use_shared_arrays;
    int switch_flag;
    double wself;
    int bzero_flag;
    std::vector<double> bzero;
};

#endif  // BISPECTRUM_HPP_
