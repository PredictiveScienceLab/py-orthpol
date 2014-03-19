================================
ORTHOGONAL POLYNOMIALS IN PYTHON
================================

Description
-----------

This package serves as a Python wrapper for the legacy Fortran code ORTHPOL.
It allows the construction of orthogonal polynomials given a weight function.


Demos
-----

I provide several demos that demonstrate how polynomials can be constructed
both from simple weight functions as well as ``scipy.stats`` random variables.
It is quite easy based on these examples to generalize to more complicated cases.
All one has to do is change the weight function or the random variable.
Here is a list of them:
+ [demos/demo1.py](demos/demo1.py): Hermite polynomials using a weight function.
+ [demos/demo2.py](demos/demo2.py): Laguerre polynomials using a weight function.
+ [demos/demo3.py](demos/demo3.py): Chebyshev polynomials using a weight function.
+ [demos/demo4.py](demos/demo4.py): Jacobi polynomials using a weight function.
+ [demos/demo5.py](demos/demo5.py): Gegenbauer polynomials using a weight function.
+ [demos/demo6.py](demos/demo6.py): Legendre polynomials using a weight function.
+ [demos/demo7.pv](demos/demo7.py): Legendre polynomials using ``scipy.stats.uniform()``.
+ [demos/demo8.pv](demos/demo8.py): Hermite polynomials using ``scipy.stats.norm()``.
+ [demos/demo9.pv](demos/demo9.py): Shifted Hermite polynomials using a non-standard ``scipy.stats.norm()``.
+ [demos/demo10.py](demos/demo10.py): Orthogonal polynomials with respect to a truncated normal.
