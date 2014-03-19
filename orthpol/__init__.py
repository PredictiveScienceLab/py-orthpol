"""Generalized Polynomial Chaos Module

Author:
    Ilias Bilionis

Date:
    8/10/2013
"""


__all__ = ['OrthogonalPolynomial', 'ProductBasis', 'QuadratureRule']


import _orthpol
from ._quadrature_rule import *
from ._orthogonal_polynomial import *
