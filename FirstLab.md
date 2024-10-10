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

---

# Exercise 2
- Spline.md
```matlab
function spline2D(knot_vectorx, knot_vectory)
    spline2Duniform(knot_vectorx, knot_vectory);
    return;
end

function spline2Duniform(knot_vectorx, knot_vectory)
    precision = 0.01;
    px = compute_p(knot_vectorx);
    nrx = compute_nr_basis_functions(knot_vectorx, px);

    py = compute_p(knot_vectory);
    nry = compute_nr_basis_functions(knot_vectory, py);

    x_begin = knot_vectorx(1);
    x_end = knot_vectorx(end);
    x = x_begin:precision:x_end;

    y_begin = knot_vectory(1);
    y_end = knot_vectory(end);
    y = y_begin:precision:y_end;
    [X, Y] = meshgrid(x, y);
    Z = zeros(size(X));

    hold on;
    for i = 1:nrx
        vx = compute_spline(knot_vectorx, px, i, X);
        for j = 1:nry
            vy = compute_spline(knot_vectory, py, j, Y);
            Z = Z + vx .* vy; 
            surf(X, Y, Z); 
        end
    end
    title('B-Spline Surface');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold off;
end

function p = compute_p(knot_vector)
    kvsize = length(knot_vector);
    i = 1;

    while (i + 1 < kvsize) && (knot_vector(i) == knot_vector(i + 1))
        i = i + 1;
    end

    p = i - 1;
    return;
end

function nr = compute_nr_basis_functions(knot_vector, p)
    nr = length(knot_vector) - p - 1; 
end

function y = compute_spline(knot_vector, p, nr, x)
    a = knot_vector(nr);
    b = knot_vector(nr + p);
    c = knot_vector(nr + 1);
    d = knot_vector(nr + p + 1);

    if (p == 0)
        y = (x >= a) & (x <= b);
        return;
    end

    lp = compute_spline(knot_vector, p - 1, nr, x);
    rp = compute_spline(knot_vector, p - 1, nr + 1, x);

    if (a == b)
        y1 = (x >= a) & (x <= b);
    else
        y1 = (x - a) / (b - a) .* (x >= a & x <= b);
    end

    if (c == d)
        y2 = (x > c) & (x <= d);
    else
        y2 = (x - c) / (d - c) .* (x > c & x <= d);
    end

    y = lp .* y1 + rp .* y2;
end
```

### Exercise 2, Task 1
- Prepare a vector of knots that generates a basis equivalent to the Lagrange basis of
polynomials of the first order on four elements [0,1], [1,2], [2,3] and [3,4] along the x
axis. Prepare a vector of knots that generates a basis equivalent to the quadratic B-
splines basis of polynomials with C1 continuity on three elements [0,1], [1,2], and
[2,3] along the y axis. Please run Octave code spline2D.m and draw base functions.
Please return the knot vectors and the 3D plot.

- main.md
```matlab

knot_vectorx = [0, 1, 2, 3, 4]; 
knot_vectory = [0, 1, 2]; 


spline2D(knot_vectorx, knot_vectory);

disp('Lagrange basis (X-axis):');
disp(knot_vectorx);
disp('B-spline basis (Y-axis):');
disp(knot_vectory);
```

<details>
<summary><b> Result 2.1 </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/dc6a4b61-0a6d-4df8-9dee-38323e4de99d">
</details>

<details>
<summary><b> Result: Command Window </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/03f2fc99-345b-499b-b711-7379e0936e41">
</details>

### Exercise 2, Task 2
- Prepare a vector of knots that generates a basis equivalent to the Lagrange basis of
polynomials of the second order on two elements [0,1], [1,2] along the x axis.
Prepare a vector of knots that generates a basis equivalent to the Lagrange basis of
polynomials of the second order on three elements [0,0.1], [0.1,0.9], and [0.9,1]
along the y axis. Please run Octave code spline2D.m and draw base functions. Please
return the knot vectors and the 3D plot

- spline2D.md
```matlab
function spline2D(knot_vectorx, knot_vectory)
    spline2Duniform(knot_vectorx, knot_vectory);
end

function spline2Duniform(knot_vectorx, knot_vectory)

    precision = 0.01;

    px = compute_p(knot_vectorx);
    nrx = compute_nr_basis_functions(knot_vectorx, px);

    py = compute_p(knot_vectory);
    nry = compute_nr_basis_functions(knot_vectory, py);

    x_begin = knot_vectorx(1);
    x_end = knot_vectorx(end);
    x = x_begin:precision:x_end;

    y_begin = knot_vectory(1);
    y_end = knot_vectory(end);
    y = y_begin:precision:y_end;
    [X, Y] = meshgrid(x, y);
    Z = zeros(size(X));

    hold on;
    for i = 1:nrx
        vx = compute_spline(knot_vectorx, px, i, X);
        for j = 1:nry
            vy = compute_spline(knot_vectory, py, j, Y);
            Z = Z + vx .* vy;
            surf(X, Y, Z); 
        end
    end
    title('B-Spline Surface');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold off;
end

function t = check_sanity(knot_vector, p)
    % Check the sanity of the knot vector
    initial = knot_vector(1);
    kvsize = size(knot_vector, 2);
    t = true;
    counter = 1;

    for i = 1:p+1
        if (initial ~= knot_vector(i))
            t = false;
            return;
        end
    end

    for i = p+2:kvsize-p-1
        if (initial == knot_vector(i))
            counter = counter + 1;
            if (counter > p)
                t = false;  
                return;
            end
        else
            initial = knot_vector(i);
            counter = 1;
        end
    end

    initial = knot_vector(kvsize);

    for i = kvsize-p:kvsize
        if (initial ~= knot_vector(i))
            t = false;
            return;
        end
    end

    for i = 1:kvsize-1
        if (knot_vector(i) > knot_vector(i+1))
            t = false;
        end
    end

    return;
end

function nr = compute_nr_basis_functions(knot_vector, p)
    nr = size(knot_vector, 2) - p - 1;
end

function p = compute_p(knot_vector)
    initial = knot_vector(1);
    kvsize = size(knot_vector, 2);
    i = 1;

    while (i+1 < kvsize) && (initial == knot_vector(i+1))
        i = i + 1;
    end
    
    p = i - 1;    
    return;
end

function y = compute_spline(knot_vector, p, nr, x)
    fC = @(x, a, b) (x - a) / (b - a);
    fD = @(x, c, d) (1 - x) / (d - c);

    a = knot_vector(nr);
    b = knot_vector(nr + p);
    c = knot_vector(nr + 1);
    d = knot_vector(nr + p + 1);

    if (p == 0)
        y = zeros(size(x));
        y(a <= x & x <= d) = 1;
        return;
    end

    lp = compute_spline(knot_vector, p - 1, nr, x);
    rp = compute_spline(knot_vector, p - 1, nr + 1, x);
    
    if (a == b)
        y1 = zeros(size(x));
        y1(a <= x & x <= b) = 1;
    else
        y1 = zeros(size(x));
        y1(a <= x & x <= b) = fC(x(a <= x & x <= b), a, b);
    end
    
    if (c == d)
        y2 = zeros(size(x));
        y2(c < x & x <= d) = 1;
    else
        y2 = zeros(size(x));
        y2(c < x & x <= d) = fD(x(c < x & x <= d), c, d);
    end

    y = lp .* y1 + rp .* y2;
    return;
end
```

- main.md
```matlab
knot_vectorx = [0, 0, 1, 2, 2];
knot_vectory = [0, 0, 0.1, 0.9, 1, 1]; 

spline2D(knot_vectorx, knot_vectory);
disp('Lagrange basis (X-axis):');
disp(knot_vectorx);
disp('Lagrange basis (Y-axis):');
disp(knot_vectory);
```

<details>
<summary><b> Result 2.2 </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/f019c6e6-5e6f-45e3-ba33-de4f9b564dfe">
</details>

<details>
<summary><b>  Command Window </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/d7950ec3-8064-4154-9c53-7f1515463d97">
</details>

### Exercise 2, Task 3
- Prepare a vector of knots that generates quadratic B-splines basis C1 on two
elements [0,1], [1,2] along the x axis. Prepare a vector of knots that generates cubic
B-splines C2 on three elements [0,1], [1,2], and [1,3] along the y axis. Please run
Octave code spline2D.m and draw base functions. Please return the knot vectors and
the 3D plot.

-spline2D.m

```matlab
function spline2D(knot_vectorx, knot_vectory)
    spline2Duniform(knot_vectorx, knot_vectory);
end

function spline2Duniform(knot_vectorx, knot_vectory)
    precision = 0.01;


    px = compute_p(knot_vectorx);
    nrx = compute_nr_basis_functions(knot_vectorx, px);

    py = compute_p(knot_vectory);
    nry = compute_nr_basis_functions(knot_vectory, py);

    x_begin = knot_vectorx(1);
    x_end = knot_vectorx(end);
    x = x_begin:precision:x_end;

    y_begin = knot_vectory(1);
    y_end = knot_vectory(end);
    y = y_begin:precision:y_end;

    [X, Y] = meshgrid(x, y);
    Z = zeros(size(X));

    hold on;
    for i = 1:nrx
        vx = compute_spline(knot_vectorx, px, i, X);
        for j = 1:nry
            vy = compute_spline(knot_vectory, py, j, Y);
            Z = Z + vx .* vy; 
            surf(X, Y, Z);
        end
    end
    title('B-Spline Surface');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold off;
end

function t = check_sanity(knot_vector, p)
    initial = knot_vector(1);
    kvsize = size(knot_vector, 2);
    t = true;
    counter = 1;

    for i = 1:p+1
        if (initial ~= knot_vector(i))
            t = false;
            return;
        end
    end

    for i = p+2:kvsize-p-1
        if (initial == knot_vector(i))
            counter = counter + 1;
            if (counter > p)
                t = false;  
                return;
            end
        else
            initial = knot_vector(i);
            counter = 1;
        end
    end

    initial = knot_vector(kvsize);

    for i = kvsize-p:kvsize
        if (initial ~= knot_vector(i))
            t = false;
            return;
        end
    end

    for i = 1:kvsize-1
        if (knot_vector(i) > knot_vector(i+1))
            t = false;
        end
    end

    return;
end

function nr = compute_nr_basis_functions(knot_vector, p)
    nr = size(knot_vector, 2) - p - 1;
end

function p = compute_p(knot_vector)
    initial = knot_vector(1);
    kvsize = size(knot_vector, 2);
    i = 1;

    while (i+1 < kvsize) && (initial == knot_vector(i+1))
        i = i + 1;
    end
    
    p = i - 1;    
    return;
end

function y = compute_spline(knot_vector, p, nr, x)
    fC = @(x, a, b) (x - a) / (b - a);
    fD = @(x, c, d) (1 - x) / (d - c);

    a = knot_vector(nr);
    b = knot_vector(nr + p);
    c = knot_vector(nr + 1);
    d = knot_vector(nr + p + 1);

    if (p == 0)
        y = zeros(size(x));
        y(a <= x & x <= d) = 1;
        return;
    end

    lp = compute_spline(knot_vector, p - 1, nr, x);
    rp = compute_spline(knot_vector, p - 1, nr + 1, x);
    
    if (a == b)
        y1 = zeros(size(x));
        y1(a <= x & x <= b) = 1;
    else
        y1 = zeros(size(x));
        y1(a <= x & x <= b) = fC(x(a <= x & x <= b), a, b);
    end
    
    if (c == d)
        y2 = zeros(size(x));
        y2(c < x & x <= d) = 1;
    else
        y2 = zeros(size(x));
        y2(c < x & x <= d) = fD(x(c < x & x <= d), c, d);
    end

    y = lp .* y1 + rp .* y2;
    return;
end
```

- main.m
  
```matlab
knot_vectorx = [0, 0, 0, 1, 1, 2, 2, 2];
knot_vectory = [0, 0, 0, 1, 1, 2, 1, 3, 3, 3]; 
spline2D(knot_vectorx, knot_vectory);

disp('Quadratic B-splines (X-axis):');
disp(knot_vectorx);
disp('Cubic B-splines (Y-axis):');
disp(knot_vectory);
```


<details>
<summary><b>  Result: 2.3 </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/794095b2-ace0-4135-9507-aa9a38bcb143">
</details>

<details>
<summary><b>  Command Window </b></summary>
<img width="745" alt="image" src="https://github.com/user-attachments/assets/99ab224d-fdf5-4aa3-94f1-db656ed4c0c7">
</details>



