#pragma once

#include <vector>
#include <matrix.hpp>

/// @brief Gauss elimination algorithm
/// @param L - matrix
/// @param P - vector with values
/// @param n_limit (default -1) - count of rows to eliminate
/// @param debug (default false) - if to print elimination completion in %
void gauss_elimination(matrix_t& L, vector_t& P, std::size_t n_limit = -1, bool debug = false);