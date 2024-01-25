#pragma once

#include <vector>
#include <matrix.hpp>

/// @brief Solves equations in triangle form
/// @param L - matrix of coefficients
/// @param P - vector of values
/// @param n_limit (default -1) - count of rows to solve
/// @param debug (default false) - if to print solve completion in %
/// @return vector of solved values
vector_t solve_triangle(matrix_t& L, vector_t& P, std::size_t n_limit = -1, bool debug = false);