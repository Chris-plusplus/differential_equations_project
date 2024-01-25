#include <solve_triangle.hpp>
#include <iostream>

vector_t solve_triangle(matrix_t& L, vector_t& P, std::size_t n_limit, bool debug) {
	std::size_t n = P.size();
	if (n_limit != (std::size_t)-1) {
		n = n_limit;
	}
	vector_t result(n);

	std::size_t last_percent = 0;
	if (debug) std::cout << "0%";
	std::size_t debug_row = 0;

	// subtracts rows upwards, starting from last row
	for (std::size_t row = n - 1; row != (std::size_t)-1; --row, ++debug_row) {
		result[row] = P[row];
		for (std::size_t row2 = row + 1; row2 != n; ++row2) {
			result[row] -= P[row2] * L[row][row2];
		}
		result[row] /= L[row][row];
		P[row] = result[row];

		if (debug and std::size_t((double)debug_row / (double)n * 100.0) != last_percent) {
			last_percent = std::size_t((double)debug_row / (double)n * 100.0);
			std::cout << '\r' << last_percent << "%";
		}
	}
	if (debug) std::cout << "\r100%\n";

	for (std::size_t i = n; i != n; ++i) {
		result[i] = P[i];
	}

	return result;
}