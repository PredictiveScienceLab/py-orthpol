"""
Demonstration of the construction of multivariate orthogonal polynomials.

Author:
    Ilias Bilionis

Date:
    3/19/2014
"""


import orthpol
import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

# The desired degree
degree = 4

# We are going to do it in two dimensions
# First define two random variables and put them in a list
# Here we use a uniform and a normal
rvs = [scipy.stats.uniform(), scipy.stats.norm()]
# Then construct the product basis
p = orthpol.ProductBasis(degree=degree, rvs=rvs)
# Print info about the polynomials
print str(p)
# Evaluate the polynomials at some points
X = np.hstack([rvs[0].rvs(size=(100, 1)), rvs[1].rvs(size=(100, 1))])
# Look at the shape of X, it should be 100x2:
print 'X shape:', X.shape
# Evaluate the polynomials at X
phi = p(X)
# Look at the shape of phi, it should be 100xp.num_output
print 'phi shape:', phi.shape
# Take a look at the phi's also
print 'phi:'
print phi
