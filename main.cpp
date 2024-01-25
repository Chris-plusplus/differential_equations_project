#include <gauss_legendre.hpp>
#include <legendre_polynomial.hpp>
#include <iostream>
#include <iomanip>
#include <gauss_elimination.hpp>
#include <solve_triangle.hpp>
#include <numbers>
#include <equation.hpp>
#include <fstream>

int main(int argc, char** argv) {
	using equation::n;

	if (argc < 3) {
		if (argc < 2) {
			std::cout << "Input finite elements count: ";
			std::cin >> n;
		}
		std::cout << "Input integration points: ";
		std::cin >> equation::integral_n;
	}
	if (argc > 1) {
		n = std::atoi(argv[1]);
		if (argc > 2) {
			equation::integral_n = std::atoi(argv[2]);
		}
	}

	std::size_t last_percent = 0;

	std::cout << "creating matrix...\n";
	matrix_t B_matrix(n + 1, vector_t(n + 1));
	vector_t L_vector(n + 1);
	std::cout << "created matrix\n";

	std::cout << "writing to matrix...\n";
	std::cout << "0%";
	for (std::size_t i = 0; i != n; ++i) {
		B_matrix[i][i] = equation::B(i, i);
		L_vector[i] = equation::L(i) - equation::B_shift(i);

		if (i + 1 != n) {
			B_matrix[i][i + 1] = equation::B(i, i + 1);
			B_matrix[i + 1][i] = B_matrix[i][i + 1];
		}

		if (std::size_t((double)i / (double)n * 100.0) != last_percent) {
			last_percent = std::size_t((double)i / (double)n * 100.0);
			std::cout << '\r' << last_percent << "%";
		}
	}
	
	/*
	for (size_t y = 0; y != n + 1; ++y) {
		for (size_t x = 0; x != n + 1; ++x) {
			std::cout << std::fixed << std::setprecision(3) << B_matrix[y][x] << (x == n ? '\n' : '\t');
		}
	}
	*/

	// Dirichlet boundary condition
	B_matrix[n][n] = 1.0;
	L_vector[n] = equation::ub;

	std::cout << "\r100%\n";
	std::cout << "written to matrix\n";

	std::cout << "eliminating with gauss method...\n";
	gauss_elimination(B_matrix, L_vector, n, true);
	std::cout << "eliminated with gauss method\n";

	std::cout << "solving...\n";
	auto result = solve_triangle(B_matrix, L_vector, n, true);
	std::cout << "solved\n";

	// desired function
	auto u = [&](double x) -> double {
		double sum = equation::shift(x);
		for (std::size_t i = 0; i != result.size(); ++i) {
			sum += equation::e(i, x) * result[i];
		}

		return sum;
	};

	std::ofstream file{"result.csv"};
	file << "x\tu(x)\n";
	std::cout << "saving results...\n";
	for (std::size_t i = 0; i != n + 1; ++i) {
		// WARNING! results are saved with dot as decimal separator
		file << equation::xi(i) << '\t' << u(equation::xi(i)) << '\n';
	}
	std::cout << "saved results\n";
}