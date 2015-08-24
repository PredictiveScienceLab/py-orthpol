"""Unit-tests for orthpol

Author:
    Ilias Bilionis

Date:
    8/10/2013

"""


import unittest
import numpy as np
import math
import orthpol
import scipy.stats as st
import matplotlib.pyplot as plt
import orthpol


class GpcTest(unittest.TestCase):

    def test_normalization(self):
        rvs = [st.uniform(loc = -1., scale = 2.)]*2
        p = orthpol.ProductBasis(rvs, degree = 5)
        U = st.uniform.rvs(loc = -1, scale = 2., size = (10000,2))
        P = p(U)
        print np.sum(P**2, axis = 0) / 10000

if __name__ == '__main__':
    unittest.main()
