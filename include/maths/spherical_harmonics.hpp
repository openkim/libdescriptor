#ifndef SPHERICAL_HARMONICS_HPP
#define SPHERICAL_HARMONICS_HPP
//From Spherical Harmonic Lighting: The Gritty Details
// No longer the case: -and google's spherical harmonic library: https://github.com/google/spherical-harmonics-
// Based on http://phys.uri.edu/nigh/NumRec/bookfpdf/f6-8.pdf Numerical Recipes in F77

//Convention, largely follows the chemist's notation, simply putting larger numbers first
// first l, then m, then theta, then phi
// theta is the polar angle, phi is the azimuthal angle
// theta is in [0, pi], phi is in [0, 2pi], m is in [-l, l]
// Condon-Shortley phase: (-1)^m

#include <iostream>
#define _USE_MATH_DEFINES

#include <cmath>
#include <stdexcept>
#include <array>
#include <vector>

#define HARD_CODED_SH_COEFF 4
#define FACT_CACHE_SIZE 16


// factorials and double facotrials ************************************************

// The vast majority of SH evaluations will hit these precomputed values.
double factorial(int x);

double factorial(double x);

double double_factorial(int x) ;
// End factorials and double factorials ********************************************

// Misc math functions *************************************************************
// clamp the first argument to be greater than or equal to the second
// and less than or equal to the third.
double clamp(double val, double min, double max) ;

// Return true if the first value is within epsilon of the second value.
bool approx_equal(double actual, double expected) ;

// Return floating mod x % m.
double fmod(double x, double m) ;

std::array<double,3> sph2cart(double phi, double theta) ;

std::array<double,3> cart2sph(const std::array<double,3> &r) ;


// *********************************************************************************

// Condon-Shortley phase: (-1)^m
double condon_shortley(int m) ;

// Evaluate Associated Legendre Polynomial P(l,m) at x
// Evaluate the associated Legendre polynomial of degree @l and order @m at
// coordinate @x. The inputs must satisfy:
// 1. l >= 0
// 2. 0 <= m <= l
// 3. -1 <= x <= 1
// See http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
//
double Plm(int l, int m, double theta) ;

double kron(int i, int j) ;

// Normalization Factor for SH Function
double K(int l, int m) ;

// Evaluate SH Function at (theta,phi) **********************************************
std::array<double,2> Ylmi(int l, int m, double phi, double theta) ;

double Ylm(int l, int m, double phi, double theta) ;

std::vector<double> Ylm_all_m(int l, double phi, double theta) ;

std::vector<double> Ylmi_all_m(int l, double phi, double theta) ;

std::vector<double> Ylm_all_m_from_r(int l, const std::array<double,3> &r) ;

std::vector<double> Ylmi_all_m_from_r(int l, const std::array<double,3> &r) ;

// pointer based version of above two functions
void Ylm_all_m(int l, double phi, double theta, double * sh) ;

void Ylmi_all_m(int l, double phi, double theta, double * sh) ;

void Ylm_all_m_from_r(int l,  double * r, double * sh) ;

void Ylmi_all_m_from_r(int l, double * r, double * sh);

void Ylmi_all_l_from_r(int l_max, double * r, double * sh);

void Ylmi_all_l_from_r(int l_max, double * r, double * sh_real, double * sh_imag);

void Ylm_all_l_from_r(int l_max, double * r, double * sh);


// All analytical SH functions *****************************************************
std::array<double, 2> Y_0_0(double theta, double phi);
std::array<double, 2> Y_1_neg1(double theta, double phi);
std::array<double, 2> Y_1_0(double theta, double phi);
std::array<double, 2> Y_1_1(double theta, double phi);
std::array<double, 2> Y_2_neg2(double theta, double phi);
std::array<double, 2> Y_2_neg1(double theta, double phi);
std::array<double, 2> Y_2_0(double theta, double phi);
std::array<double, 2> Y_2_1(double theta, double phi);
std::array<double, 2> Y_2_2(double theta, double phi);
std::array<double, 2> Y_3_neg3(double theta, double phi);
std::array<double, 2> Y_3_neg2(double theta, double phi);
std::array<double, 2> Y_3_neg1(double theta, double phi);
std::array<double, 2> Y_3_0(double theta, double phi);
std::array<double, 2> Y_3_1(double theta, double phi);
std::array<double, 2> Y_3_2(double theta, double phi);
std::array<double, 2> Y_3_3(double theta, double phi);
std::array<double, 2> Y_4_neg4(double theta, double phi);
std::array<double, 2> Y_4_neg3(double theta, double phi);
std::array<double, 2> Y_4_neg2(double theta, double phi);
std::array<double, 2> Y_4_neg1(double theta, double phi);
std::array<double, 2> Y_4_0(double theta, double phi);
std::array<double, 2> Y_4_1(double theta, double phi);
std::array<double, 2> Y_4_2(double theta, double phi);
std::array<double, 2> Y_4_3(double theta, double phi);
std::array<double, 2> Y_4_4(double theta, double phi);
std::array<double, 2> Y_5_neg5(double theta, double phi);
std::array<double, 2> Y_5_neg4(double theta, double phi);
std::array<double, 2> Y_5_neg3(double theta, double phi);
std::array<double, 2> Y_5_neg2(double theta, double phi);
std::array<double, 2> Y_5_neg1(double theta, double phi);
std::array<double, 2> Y_5_0(double theta, double phi);
std::array<double, 2> Y_5_1(double theta, double phi);
std::array<double, 2> Y_5_2(double theta, double phi);
std::array<double, 2> Y_5_3(double theta, double phi);
std::array<double, 2> Y_5_4(double theta, double phi);
std::array<double, 2> Y_5_5(double theta, double phi);
std::array<double, 2> Y_6_neg6(double theta, double phi);
std::array<double, 2> Y_6_neg5(double theta, double phi);
std::array<double, 2> Y_6_neg4(double theta, double phi);
std::array<double, 2> Y_6_neg3(double theta, double phi);
std::array<double, 2> Y_6_neg2(double theta, double phi);
std::array<double, 2> Y_6_neg1(double theta, double phi);
std::array<double, 2> Y_6_0(double theta, double phi);
std::array<double, 2> Y_6_1(double theta, double phi);
std::array<double, 2> Y_6_2(double theta, double phi);
std::array<double, 2> Y_6_3(double theta, double phi);
std::array<double, 2> Y_6_4(double theta, double phi);
std::array<double, 2> Y_6_5(double theta, double phi);
std::array<double, 2> Y_6_6(double theta, double phi);
std::array<double, 2> Y_7_neg7(double theta, double phi);
std::array<double, 2> Y_7_neg6(double theta, double phi);
std::array<double, 2> Y_7_neg5(double theta, double phi);
std::array<double, 2> Y_7_neg4(double theta, double phi);
std::array<double, 2> Y_7_neg3(double theta, double phi);
std::array<double, 2> Y_7_neg2(double theta, double phi);
std::array<double, 2> Y_7_neg1(double theta, double phi);
std::array<double, 2> Y_7_0(double theta, double phi);
std::array<double, 2> Y_7_1(double theta, double phi);
std::array<double, 2> Y_7_2(double theta, double phi);
std::array<double, 2> Y_7_3(double theta, double phi);
std::array<double, 2> Y_7_4(double theta, double phi);
std::array<double, 2> Y_7_5(double theta, double phi);
std::array<double, 2> Y_7_6(double theta, double phi);
std::array<double, 2> Y_7_7(double theta, double phi);
std::array<double, 2> Y_8_neg8(double theta, double phi);
std::array<double, 2> Y_8_neg7(double theta, double phi);
std::array<double, 2> Y_8_neg6(double theta, double phi);
std::array<double, 2> Y_8_neg5(double theta, double phi);
std::array<double, 2> Y_8_neg4(double theta, double phi);
std::array<double, 2> Y_8_neg3(double theta, double phi);
std::array<double, 2> Y_8_neg2(double theta, double phi);
std::array<double, 2> Y_8_neg1(double theta, double phi);
std::array<double, 2> Y_8_0(double theta, double phi);
std::array<double, 2> Y_8_1(double theta, double phi);
std::array<double, 2> Y_8_2(double theta, double phi);
std::array<double, 2> Y_8_3(double theta, double phi);
std::array<double, 2> Y_8_4(double theta, double phi);
std::array<double, 2> Y_8_5(double theta, double phi);
std::array<double, 2> Y_8_6(double theta, double phi);
std::array<double, 2> Y_8_7(double theta, double phi);
std::array<double, 2> Y_8_8(double theta, double phi);
std::array<double, 2> Y_9_neg9(double theta, double phi);
std::array<double, 2> Y_9_neg8(double theta, double phi);
std::array<double, 2> Y_9_neg7(double theta, double phi);
std::array<double, 2> Y_9_neg6(double theta, double phi);
std::array<double, 2> Y_9_neg5(double theta, double phi);
std::array<double, 2> Y_9_neg4(double theta, double phi);
std::array<double, 2> Y_9_neg3(double theta, double phi);
std::array<double, 2> Y_9_neg2(double theta, double phi);
std::array<double, 2> Y_9_neg1(double theta, double phi);
std::array<double, 2> Y_9_0(double theta, double phi);
std::array<double, 2> Y_9_1(double theta, double phi);
std::array<double, 2> Y_9_2(double theta, double phi);
std::array<double, 2> Y_9_3(double theta, double phi);
std::array<double, 2> Y_9_4(double theta, double phi);
std::array<double, 2> Y_9_5(double theta, double phi);
std::array<double, 2> Y_9_6(double theta, double phi);
std::array<double, 2> Y_9_7(double theta, double phi);
std::array<double, 2> Y_9_8(double theta, double phi);
std::array<double, 2> Y_9_9(double theta, double phi);
std::array<double, 2> Y_10_neg10(double theta, double phi);
std::array<double, 2> Y_10_neg9(double theta, double phi);
std::array<double, 2> Y_10_neg8(double theta, double phi);
std::array<double, 2> Y_10_neg7(double theta, double phi);
std::array<double, 2> Y_10_neg6(double theta, double phi);
std::array<double, 2> Y_10_neg5(double theta, double phi);
std::array<double, 2> Y_10_neg4(double theta, double phi);
std::array<double, 2> Y_10_neg3(double theta, double phi);
std::array<double, 2> Y_10_neg2(double theta, double phi);
std::array<double, 2> Y_10_neg1(double theta, double phi);
std::array<double, 2> Y_10_0(double theta, double phi);
std::array<double, 2> Y_10_1(double theta, double phi);
std::array<double, 2> Y_10_2(double theta, double phi);
std::array<double, 2> Y_10_3(double theta, double phi);
std::array<double, 2> Y_10_4(double theta, double phi);
std::array<double, 2> Y_10_5(double theta, double phi);
std::array<double, 2> Y_10_6(double theta, double phi);
std::array<double, 2> Y_10_7(double theta, double phi);
std::array<double, 2> Y_10_8(double theta, double phi);
std::array<double, 2> Y_10_9(double theta, double phi);
std::array<double, 2> Y_10_10(double theta, double phi);
std::array<double, 2> Ylmi_analytical(int l,int m, double theta, double phi);

#endif
