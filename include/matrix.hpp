#pragma once

#include <vector>
#include <unordered_map>
#include <ranges>
#include <algorithm>

/// @brief Vector of doubles
using vector_t = std::vector<double>;
/// @brief Square matrix of doubles
/// @brief WARNING! martix is not sparse, might use a lot of memory (ex. n = 30000 uses >4GB of RAM)
using matrix_t = std::vector<vector_t>;

namespace mx {
/// @brief Adds rows of matrix, along with associated value from vector
/// @param L - matrix
/// @param P - vector
/// @param to_add_to - row index of row to add to
/// @param added - row index of row to add
/// @param scalar (default = 1.0) - value to multiply added row with
void add_rows(matrix_t& L, vector_t& P, std::size_t to_add_to, std::size_t added, double scalar = 1.0);
/// @brief Subtracts rows of matrix, along with associated value from vector
/// @param L - matrix
/// @param P - vector
/// @param to_add_to - row index of row to subtract from
/// @param subtracted - row index of row to subtract
/// @param scalar (default = 1.0) - value to multiply subtracted row with
void subtract_rows(matrix_t& L, vector_t& P, std::size_t to_subtract_from, std::size_t subtracted, double scalar = 1.0);
}