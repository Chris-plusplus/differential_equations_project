#include <gauss_legendre.hpp>
#include <legendre_polynomial.hpp>
#include <algorithm>
#include <numbers>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <cmath>

legendre_polynomial::legendre_polynomial(std::size_t n):
	_roots(n),
	_weights(n) {
	compute();
}

const legendre_polynomial& legendre_polynomial::get(std::size_t n) {
	// searches for cached
	auto found = _computed_polys.find(n);
	if (found != _computed_polys.end()) {
		// found cached polynomial
		return found->second;
	}
	else {
		// not found cached polynomial

		// searches for file-cached polynomial
		if (std::filesystem::exists("cached_legendre_polynomials/" + std::to_string(n))) {
			// found file-cached polynomial

			legendre_polynomial lpoly{};
			lpoly._roots = std::vector<double>(n);
			lpoly._weights = std::vector<double>(n);

			// file cache is read as raw double array
			// 100% safe xd
			std::ifstream file{"cached_legendre_polynomials/" + std::to_string(n), std::ios::binary | std::ios::in};
			file.read((char*)lpoly._roots.data(), n * sizeof(double));
			file.read((char*)lpoly._weights.data(), n * sizeof(double));

			_computed_polys.insert_or_assign(n, lpoly);

			return get(n);
		}
		else {
			// not found file-cached polynomial

			legendre_polynomial lpoly{n};
			if (not std::filesystem::exists("cached_legendre_polynomials")) {
				std::filesystem::create_directory("cached_legendre_polynomials");
			}
			
			// file cache is saved as raw double array
			// 100% safe xd
			std::ofstream file{"cached_legendre_polynomials/" + std::to_string(n), std::ios::binary | std::ios::out};
			file.write((char*)lpoly._roots.data(), n * sizeof(double));
			file.write((char*)lpoly._weights.data(), n * sizeof(double));

			_computed_polys.insert_or_assign(n, lpoly);

			return get(n);
		}
	}
}

std::size_t legendre_polynomial::n() const {
	return _roots.size();
}
const std::vector<double>& legendre_polynomial::weights() const {
	return _weights;
}
const std::vector<double>& legendre_polynomial::roots() const {
	return _roots;
}

legendre_polynomial::value_derivative_pair legendre_polynomial::compute_vd(double x) const {
	// value of 1st (not 0th) legendre polynomial
	vd_pair result{x, 0.0};

	double value_old = 1.0;
	for (std::size_t i = 2; i <= n(); ++i) {
		// recursive formula for polynomial value
		double ith_value = ((2 * i - 1) * x * result.value - (i - 1) * value_old) / i;

		value_old = result.value;
		result.value = ith_value;
	}

	result.derivative = n() * 1/(x*x-1) * (x * result.value - value_old);

	return result;
}

void legendre_polynomial::compute() {
	for (std::size_t i = 0; i != n(); ++i) {
		// root approximation formula available on German wikipedia
		double root = cos(std::numbers::pi_v<double> *((double)i - 0.25) / ((double)n() + 0.5));
		
		auto vd = compute_vd(root);
		double fdf_ratio;
		do { // computing roots using Newton's method
			fdf_ratio = vd.value / vd.derivative;
			root -= fdf_ratio;
			vd = compute_vd(root);
		}
		while (std::fabs(fdf_ratio) > epsilon);

		_roots[n() - i - 1] = root;
		_weights[n() - i - 1] = 2.0 / ((1 - root * root) * vd.derivative * vd.derivative);
	}
}