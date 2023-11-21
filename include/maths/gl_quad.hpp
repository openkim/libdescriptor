#ifndef GL_QUAD_HPP
#define GL_QUAD_HPP

#include <vector>

// gauss-Lebedev quadrature for size 100, same as taken by dscribe
// TODO: higher accuracy quadrature?

std::vector<double> get_gl_weights();
std::vector<double> get_orig_gl_grid();

std::vector<double> get_gl_grid(double cutoff);
#endif //GL_QUAD_HPP