# Elastic Deformation FEM project

Computes Elastic Deformation function using **[Finite Element Method](https://en.wikipedia.org/wiki/Finite_element_method)**.

Written in `C++20` as final project for Differential and Difference Equations course.

Graded 50/50 pts

## Problem
Find u(x) satisfying:

<img src="img/problem_desc.png" width="600">

For 50 pts, it was required for all calculations (beside formula derivation) to be computed numerically.

Integrals had to be computed using **[Gauss-Legendre quadrature](https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature)** using a minimum of 2 points.

## Solution

### Derivation of equation and matrix
<img src="1.jpg" width="750">
<img src="2.jpg" width="750">

### Computations

**Gauss-Legendre quadrature** is computed using formula:

<img src="img/integral.png" width="750">

where:

<img src="img/root_in_formula.png" width="25"> - i-th root of n-th Legendre polynomial

<img src="img/weight_in_formula.png" width="25"> - i-th weight, computed with formula:

<img src="img/weight.png" width="300">

My solution computes its own roots of **[Legendre Polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials)** using **[Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method)**

<img src="img/newton.png" width="300">

where initial approximation formula is:

(found on German Wikipedia of **[Legendre Polynomials](https://de.wikipedia.org/wiki/Legendre-Polynom#Nullstellen)**)

<img src="img/root.png" width="450">

*f*(x) and *f*'(x) were computed using recursive formulas:

<img src="img/legendre_poly_recursive_formula.png" width="450">

<img src="img/legendre_poly_derivative_recursive_formula.png" width="450">

Once computed, n-th Legendre polynomial is saved to file `./cached_legendre_polynomials/<n>`, to preserve time

(~16 minutes to compute roots of 200000th polynomial)

Systems of equations shown in second picture are then solved using **[Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination)** and are saved to file `result.csv`

## Building
### Clone repo
```bash
git clone https://github.com/Chris-plusplus/differential_equations_project.git
```
### Build using CMake
```bash
cd differential_equations_project
mkdir build
cd build
cmake ..
cmake --build .
```

## Usage

Solution might be executed with following command:
```bash
./differential_equations_project <finite_element_count> <integral_points>
```
or (requires input of `finite_element_count` `integral_points`)
```bash
./differential_equations_project
```

where `intergral_points` is equal to index of Legendre polynomial used.

### Final solution
<img src="img/graph.png">
