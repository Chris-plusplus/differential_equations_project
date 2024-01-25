#include <matrix.hpp>

namespace mx {
void add_rows(matrix_t& L, vector_t& P, std::size_t to_add_to, std::size_t added, double scalar) {
	for (std::size_t i = 0; i != P.size(); ++i) {
		L[to_add_to][i] += scalar * L[added][i];
	}
	P[to_add_to] += scalar * P[added];
};

void subtract_rows(matrix_t& L, vector_t& P, std::size_t to_subtracted_from, std::size_t subtracted, double scalar) {
	add_rows(L, P, to_subtracted_from, subtracted, -scalar);
};
}