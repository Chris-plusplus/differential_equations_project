#pragma once

#include <legendre_polynomial.hpp>
#include <concepts>

/// @brief Gauss-Legendre integration function
/// @param fn - object satisfying std::invocable<double> (ex. function pointer or lambda)
/// @param a - integration interval begin
/// @param b - integration interval end
/// @param n (default = 1000) - points of Legendre polynomial to use
template<std::invocable<double> Fn>
double gauss_legendre_intergrate(Fn fn, double a, double b, std::size_t n = 1000) {
	double sum = 0;
	double mean = (a + b) / 2.0;
	double half_width = (b - a) / 2.0;
	auto lpoly = legendre_polynomial::get(n);

	for (std::size_t i = 0; i != n; ++i) {
		sum += lpoly.weights()[i] * fn(half_width * lpoly.roots()[i] + mean);
	}

	return half_width * sum;
}