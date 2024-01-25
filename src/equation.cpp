#include <equation.hpp>
#include <gauss_legendre.hpp>
#include <iostream>
#include <cmath>

namespace equation {
double xi(std::size_t i) {
	return (b - a) / (double)n * (double)i + a;
}

double xin(std::size_t i, std::size_t _n) {
	return (b - a) / (double)_n * (double)i + a;
}

double e(std::size_t i, double x) {
	if (i == 0) {
		if (xi(i) <= x and x <= xi(i + 1)) {
			return (xi(i + 1) - x) / (xi(i + 1) - xi(i));
		}
		else {
			return 0;
		}
	}
	else if (i == n + 1) {
		if (xi(i - 1) <= x and x <= xi(i)) {
			return (x - xi(i - 1)) / (xi(i) - xi(i - 1));
		}
		else {
			return 0;
		}
	}
	else {
		if (xi(i) <= x and x <= xi(i + 1)) {
			return (xi(i + 1) - x) / (xi(i + 1) - xi(i));
		}
		else if (xi(i - 1) <= x and x <= xi(i)) {
			return (x - xi(i - 1)) / (xi(i) - xi(i - 1));
		}
		else {
			return 0;
		}
	}
}
double de(std::size_t i, double x) {
	if (i == 0) {
		if (xi(i) <= x and x <= xi(i + 1)) {
			return -1.0 / (xi(i + 1) - xi(i));
		}
		else {
			return 0;
		}
	}
	else if (i == n + 1) {
		if (xi(i - 1) <= x and x <= xi(i)) {
			return 1.0 / (xi(i) - xi(i - 1));
		}
		else {
			return 0;
		}
	}
	else {
		if (xi(i) <= x and x <= xi(i + 1)) {
			return -1.0 / (xi(i + 1) - xi(i));
		}
		else if (xi(i - 1) <= x and x <= xi(i)) {
			return 1.0 / (xi(i) - xi(i - 1));
		}
		else {
			return 0;
		}
	}
}
double shift(double x) {
	return ub * e(n + 1, x);
}
double dshift(double x) {
	return ub * de(n + 1, x);
}

double B(std::size_t i, std::size_t j) {
	auto integrand = [&](double x) -> double {
		return de(i, x) * de(j, x) * E(x);
	};

	return gauss_legendre_intergrate(integrand, a, b, integral_n) - e(i, a) * e(j, a) * E(a) * q / p;
}

double L(std::size_t i) {
	return -E(a) * e(i, a) * r / p;
}
double B_shift(std::size_t i) {
	if (std::fabs(ub) < 1.0e-15) {
		return 0;
	}

	auto integrand = [&](double x) -> double {
		return dshift(x) * de(i, x) * E(x);
	};

	return gauss_legendre_intergrate(integrand, a, b, integral_n) - e(i, a) * shift(a) * E(a) * q / p;
}
}