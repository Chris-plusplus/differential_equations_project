#include <gauss_elimination.hpp>
#include <iostream>

void gauss_elimination(matrix_t& L, vector_t& P, std::size_t n_limit, bool debug) {
	std::size_t n = P.size();
	if (n_limit != (std::size_t)-1) {
		n = n_limit;
	}

	std::size_t last_percent = 0;
	if (debug) std::cout << "0%";
	for (std::size_t row = 0; row != n; ++row) {
		// subtracting only next row, because matrix is three-diagonal
		if (std::size_t to_subtract_from = row + 1; to_subtract_from != n) {
			double scalar = L[to_subtract_from][row] / L[row][row];
			
			mx::subtract_rows(L, P, to_subtract_from, row, scalar);
		}

		if (debug and std::size_t((double)row / (double)n * 100.0) != last_percent) {
			last_percent = std::size_t((double)row / (double)n * 100.0);
			std::cout << '\r' << last_percent << "%";
		}
	}
	if (debug) std::cout << "\r100%\n";
}