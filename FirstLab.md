# Lagrange Basis Functions

## Overview

**Lagrange Basis Functions** are fundamental polynomials used for interpolating a set of given data points. They form the foundation of **Lagrange Interpolation**, where the goal is to construct a polynomial that passes exactly through a set of defined points (or knots). 

Lagrange basis functions are used in various fields such as numerical analysis, computer-aided design (CAD), and finite element methods (FEM).

### Problem Setup

Given $$\( n \)$$ data points $$\( (x_1, y_1), (x_2, y_2)$$, $$\dots$$, $$(x_n, y_n) \)$$, the objective is to find a polynomial $$\( P(x) \)$$ of degree $$\( n-1 \)$$ that passes through these points.

The **Lagrange polynomial** \( P(x) \) is given as:

$$
P(x) = \sum_{i=1}^{n} y_i L_i(x)
$$

Where:
- $$\( L_i(x) \)$$ is the **i-th Lagrange basis function**, which is constructed such that:
  - $$\( L_i(x_j) = 1 \)$$ if $$\( i = j \)$$
  - $$\( L_i(x_j) = 0 \)$$ if $$\( i \neq j \)$$

### Lagrange Basis Function Definition

The **Lagrange basis function** $$\( L_i(x) \)$$ for the $$\( i \)-th$$ knot is defined as:

$$
L_i(x) = \prod_{\substack{1 \leq j \leq n \\ j \neq i}} \frac{x - x_j}{x_i - x_j}
$$

Where:
- $$\( x_1, x_2, \dots, x_n \)$$ are the knots or interpolation points.
- $$\( L_i(x) \) is 1 at \( x_i \)$$ and 0 at all other $$\( x_j \)$$ where $$\( j \neq i \)$$.

### Example

For 3 points $$\( x_1 = 0 \)$$, $$\( x_2 = 1 \)$$, and $$\( x_3 = 2 \)$$, the Lagrange basis functions would be:

$$
L_1(x) = \frac{(x - x_2)(x - x_3)}{(x_1 - x_2)(x_1 - x_3)} = \frac{(x - 1)(x - 2)}{(0 - 1)(0 - 2)}
$$

$$
L_2(x) = \frac{(x - x_1)(x - x_3)}{(x_2 - x_1)(x_2 - x_3)} = \frac{(x - 0)(x - 2)}{(1 - 0)(1 - 2)}
$$

$$
L_3(x) = \frac{(x - x_1)(x - x_2)}{(x_3 - x_1)(x_3 - x_2)} = \frac{(x - 0)(x - 1)}{(2 - 0)(2 - 1)}
$$

These basis functions can then be used to construct the interpolating polynomial

$$
\( P(x) \)
$$

- **Numerical Interpolation**: Used to estimate unknown values from known data points.
- **Finite Element Methods (FEM)**: Lagrange polynomials are widely used in FEM to approximate solutions to partial differential equations.
- **Computer-Aided Design (CAD)**: In CAD systems, Lagrange polynomials help in constructing curves that pass through a set of given points.


---

### Exercise 1, Task 1
- Create a vector of knots for Lagrange basis of fourth-degree polynomials on five intervals: [0,1], [1,2], [2,3], [3,4], and [4,5].

```matlab
knots = [0, 1, 2, 3, 4, 5];
syms x;
LagrangeBasis = sym(zeros(5, 5));
for i = 1:5
    Li = 1; 
    for j = 1:5
        if i ~= j
            Li = Li * (x - knots(j)) / (knots(i) - knots(j));
        end
    end
    LagrangeBasis(i, :) = Li; 
end

fplot(LagrangeBasis(1), [0, 5]);
hold on;
for i = 2:5
    fplot(LagrangeBasis(i), [0, 5]);
end
title('Lagrange Basis Functions');
xlabel('x');
ylabel('L(x)');
legend('L1', 'L2', 'L3', 'L4', 'L5');
hold off;
```

<details>
<summary><b> Result: Exercise 1, Task 1  </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/50dbcea7-e38a-433b-97e2-267dc9479f35">
</details>

---

### Exercise 1, Task 2
- Create a vector of knots for the B-spline basis of third-order with C2 continuity over the intervals: [0, 0.1], [0.1, 0.9], [0.9, 1].

```matlab
knots = [0 0 0 0 0.1 0.9 1 1 1 1]; % b-spline 

bspline_degree = 3; 
n = numel(knots) - bspline_degree - 1; 

figure;
hold on;
for i = 1:n
    bspline = spmak(knots, eye(n)); 
    fnplt(bspline, 2);
end
title('B-spline Basis Functions (Cubic, C2 Continuity)');
xlabel('x');
ylabel('B(x)');
hold off;
```

<details>
<summary><b> Result: Exercise 1, Task 2 </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/b9947f9d-39fa-4f53-b5ff-e88e3487662a">
</details>

---

### Exercise 1, Task 3
- Create a vector of knots for the B-spline basis of third-order with C1 continuity over four intervals: [0,1], [1,2], [2,3], [3,4].
  
```matlab
knots = [0 0 0 0 1 2 3 4 4 4 4]; 

bspline_degree = 3;
n = numel(knots) - bspline_degree - 1;

figure;
hold on;
for i = 1:n
    bspline = spmak(knots, eye(n));
    fnplt(bspline, 2);
end
title('B-spline Basis Functions (Cubic, C1 Continuity)');
xlabel('x');
ylabel('B(x)');
hold off;
```
<details>
<summary><b> Result: Exercise 1, Task 3 </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/981dd9bc-25a1-4742-9842-9e4bd1109014">
</details>

---

### Exercise 1, Task 4
- Create a vector of knots for the Lagrange basis of first-order polynomials on two intervals: [0,1], [1,2].
```matlab
knots = [0, 1, 2];
syms x;

LagrangeBasis = sym(zeros(2, 2));
for i = 1:2
    Li = 1;
    for j = 1:2
        if i ~= j
            Li = Li * (x - knots(j)) / (knots(i) - knots(j));
        end
    end
    LagrangeBasis(i, :) = Li;
end

fplot(LagrangeBasis(1), [0, 2]);
hold on;
fplot(LagrangeBasis(2), [0, 2]);
title('Lagrange Basis Functions (1st order)');
xlabel('x');
ylabel('L(x)');
legend('L1', 'L2');
hold off;
```

<details>
<summary><b> Result: Exercise 1, Task 4 </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/acc66709-2181-420f-a802-1e0d308689d6">
</details>

