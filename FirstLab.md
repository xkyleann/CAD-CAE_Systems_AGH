#Lagrange Basis Functions

## Overview

**Lagrange Basis Functions** are fundamental polynomials used for interpolating a set of given data points. They form the foundation of **Lagrange Interpolation**, where the goal is to construct a polynomial that passes exactly through a set of defined points (or knots). 

Lagrange basis functions are used in various fields such as numerical analysis, computer-aided design (CAD), and finite element methods (FEM).

### Problem Setup

Given \( n \) data points \( (x_1, y_1), (x_2, y_2), \dots, (x_n, y_n) \), the objective is to find a polynomial \( P(x) \) of degree \( n-1 \) that passes through these points.

The **Lagrange polynomial** \( P(x) \) is given as:

\[
P(x) = \sum_{i=1}^{n} y_i L_i(x)
\]

Where:
- \( L_i(x) \) is the **i-th Lagrange basis function**, which is constructed such that:
  - \( L_i(x_j) = 1 \) if \( i = j \)
  - \( L_i(x_j) = 0 \) if \( i \neq j \)

### Lagrange Basis Function Definition

The **Lagrange basis function** \( L_i(x) \) for the \( i \)-th knot is defined as:

\[
L_i(x) = \prod_{\substack{1 \leq j \leq n \\ j \neq i}} \frac{x - x_j}{x_i - x_j}
\]

Where:
- \( x_1, x_2, \dots, x_n \) are the knots or interpolation points.
- \( L_i(x) \) is 1 at \( x_i \) and 0 at all other \( x_j \) where \( j \neq i \).

### Example

For 3 points \( x_1 = 0 \), \( x_2 = 1 \), \( x_3 = 2 \), the Lagrange basis functions would be:

\[
L_1(x) = \frac{(x - x_2)(x - x_3)}{(x_1 - x_2)(x_1 - x_3)} = \frac{(x - 1)(x - 2)}{(0 - 1)(0 - 2)}
\]
\[
L_2(x) = \frac{(x - x_1)(x - x_3)}{(x_2 - x_1)(x_2 - x_3)} = \frac{(x - 0)(x - 2)}{(1 - 0)(1 - 2)}
\]
\[
L_3(x) = \frac{(x - x_1)(x - x_2)}{(x_3 - x_1)(x_3 - x_2)} = \frac{(x - 0)(x - 1)}{(2 - 0)(2 - 1)}
\]

These basis functions can then be used to construct the interpolating polynomial \( P(x) \).

### Applications

- **Numerical Interpolation**: Used to estimate unknown values from known data points.
- **Finite Element Methods (FEM)**: Lagrange polynomials are widely used in FEM to approximate solutions to partial differential equations.
- **Computer-Aided Design (CAD)**: In CAD systems, Lagrange polynomials help in constructing curves that pass through a set of given points.

---
### Exercise 1, Task 1
- Create a vector of knots for Lagrange basis of fourth-degree polynomials on five intervals: [0,1], [1,2], [2,3], [3,4], and [4,5].
- 
