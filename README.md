Orthogonal Polynomials in Python
================================

Description
-----------
The ``py-orthpol`` package defines the model ``orthpol`` which can be used
easily construct univariate and multivariate orthogonal polynomials in Python.
The purpose of this code is to serve as a component in Python packages that
could use orthogonal polynomials as basis functions for various tasks.
For example:
+ The polynomials can be used in least squares applications.
+ The polynomials can serve as the mean in Gaussian process regression.
+ etc.

The need to have an easy to use package that can generate polynomials orthogonal
with respect to arbitrary weight functions is motivated by applications in the
field of Uncertainty Quantification (UQ). In UQ, collections of such polynomials
are known as generalized Polynomial Chaos (gPC). My end goal is to provide a tool
that makes it **ridiculously easy** to construct these polynomials.

Where does this come from?
--------------------------

This package serves as a Python wrapper for the legacy Fortran code
[ORTHPOL](http://dl.acm.org/citation.cfm?id=174605). The original ORTHPOL
code can be found
[here](https://www.cs.purdue.edu/archives/2001/wxg/codes/ORTHPOL).
The code that computes tensors products of univariate orthogonal polynomials
is a transolation of [Stockos](http://trilinos.sandia.gov/packages/stokhos/)
C++ routines to Python.

Installation
------------

Simply clone the repository:

```
git clone https://github.com/ebilionis/py-orthpol.git
```

Go inside the directory and run:

```
python setup.py install
```

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
+ [demos/demo11.py](demos/demo11.py): 2D orthogonal polynomials using the ``ProductBasis`` class and a collection of ``scipy.stats`` random variables.


TODO
----

This is a list of things that need to be done:
+ Implement the method ``orthopol.ProductBasis.d()`` that calculates the
derivative of a product basis of polynomials with respect to x.
