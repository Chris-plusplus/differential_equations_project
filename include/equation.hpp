#pragma once
#include <type_traits>

namespace equation {
/// @brief finite element count
inline std::size_t n;
/// @brief integration iterations
inline std::size_t integral_n = 10000;

/// @brief Interval start
constexpr const double a = 0;
/// @brief Interval end
constexpr const double b = 2;

/// @brief E(x), as described on pictures 
inline double E(double x) {
	if (x <= 1.0) {
		return 3;
	}
	else {
		return 5;
	}
}

/// @brief Dirichlet boundary condition 
constexpr const double ub = 0;

/// @brief Cauchy boundary condition, as described on pictures
constexpr const double p = 1;
/// @brief Cauchy boundary condition, as described on pictures
constexpr const double q = 1;
/// @brief Cauchy boundary condition, as described on pictures
constexpr const double r = 10;

/// @brief Returns i-th x on interval [a, b]
/// @brief xi(0) = a
/// @brief xi(n) = b
double xi(std::size_t i);

/// @brief Returns value of i-th finite element for x
double e(std::size_t i, double x);
/// @brief Returns derivative value of i-th finite element for x
double de(std::size_t i, double x);
/// @brief Shift function for nonzero Dirichlet boundary condition
double shift(double x);
/// @brief Shift function derivative for nonzero Dirichlet boundary condition
double dshift(double x);

/// @brief Funcion B(u, v) (or B(w, v)), as described on pictures
/// @brief corresponds to B(ei, ej), as described on pictures
double B(std::size_t i, std::size_t j);

/// @brief Funcion L(v), as described on pictures
/// @brief corresponds to L(ei), as described on pictures
double L(std::size_t i);
/// @brief Funcion B(u_shift, v), as described on pictures
/// @brief corresponds to B(u_shift, ei), as described on pictures
double B_shift(std::size_t i);
}
