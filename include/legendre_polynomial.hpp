#pragma once

#include <vector>
#include <utility>
#include <unordered_map>

/// @brief Class representing Legendre polynomial
/// @brief All computed polynomials are cached in memory and in './cached_legendre_polynomials/<index_of_polynomial>'
class legendre_polynomial {
public:
	/// @brief Returns index of this polynomial
	std::size_t n() const;
	/// @brief Returns weigths of this polynomial
	const std::vector<double>& weights() const;
	/// @brief Returns weigths of this polynomial
	const std::vector<double>& roots() const;

	/// @brief Returns nth Legendre polynomial as legendre_polynomial object
	static const legendre_polynomial& get(std::size_t n);

	/// @brief DO NOT USE!, used by cache
	legendre_polynomial() = default;
private:
	/// @brief Computes n-th legendre polynomial roots and weigths
	legendre_polynomial(std::size_t n);

	/// @brief Values below are considered 0
	static inline constexpr double epsilon = 1.0e-15;

	/// @brief Cache
	static inline std::unordered_map<std::size_t, legendre_polynomial> _computed_polys;

	std::vector<double> _weights;
	std::vector<double> _roots;

	using vd_pair = struct value_derivative_pair {
		double value;
		double derivative;
	};

	/// @brief Computes roots and weigths
	void compute();
	/// @brief Computes value of polynomial and its derivative for x
	vd_pair compute_vd(double x) const;
};