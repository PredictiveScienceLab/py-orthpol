"""Unit-tests for best.gpc

Author:
    Ilias Bilionis

Date:
    8/10/2013

"""


import unittest
import numpy as np
import math
import best.gpc
import scipy.stats as stats
import matplotlib.pyplot as plt


class GpcTest(unittest.TestCase):

    def test_quadrature_rule(self):
        return
        # Testing [-1, 1]
        quad = best.gpc.QuadratureRule(-1, 1, ncap=10)
        print str(quad)
        # Testing [-1, 1]
        quad = best.gpc.QuadratureRule(0, 1, ncap=10)
        print str(quad)
        # Testing [-inf, inf]
        quad = best.gpc.QuadratureRule(-float('inf'), float('inf'), ncap=10)
        print str(quad)
        # Now some tests with non-identity weight functions.
        wf = lambda(x): np.exp(-x)
        # Testing [-1, 1]
        quad = best.gpc.QuadratureRule(-1, 1, wf=wf, ncap=10)
        print str(quad)
        # Testing [-1, 1]
        quad = best.gpc.QuadratureRule(0, 1, ncap=10)
        print str(quad)
        # Testing [-inf, inf]
        quad = best.gpc.QuadratureRule(-float('inf'), float('inf'), ncap=10)
        print str(quad)
        print 'We will integrate some functions'
        f = lambda(x): np.sin(x)
        quad = best.gpc.QuadratureRule(0, math.pi, ncap=50)
        print quad.integrate(f)

    def test_hermite(self):
        wf = lambda(x): 1. / math.sqrt(2. * math.pi) * np.exp(-x ** 2 / 2.)
        infty = float('inf')
        p = best.gpc.OrthogonalPolynomial(10, left=-infty, right=infty,
                                          wf=wf)
        print p._evaluate_square_norms()
        print str(p)
        x = np.linspace(-2., 2., 100)
        phi = p(x)
        dphi = p.d(x)

    def test_laguerre(self):
        wf = lambda(x): np.exp(-x)
        infty = float('inf')
        p = best.gpc.OrthogonalPolynomial(10, left=0, right=infty,
                                          wf=wf)
        print p._evaluate_square_norms()
        print str(p)
        import scipy.stats
        rv = scipy.stats.expon()
        p = best.gpc.OrthogonalPolynomial(10, left=0, right=infty,
                                          wf=rv.pdf)
        x = np.linspace(0., 5., 100)
        phi = p(x)
        dphi = p.d(x)

    def test_orthogonal_polynomial(self):
        import time
        p = best.gpc.OrthogonalPolynomial(3, left=-1, right=1, wf=lambda(x): 0.5, ncap=50)
        print str(p)
        print p._evaluate_square_norms()
        x = np.linspace(-1, 1, 100)
        start_time = time.time()
        phi = p(x)
        end_time = time.time()
        print("Elapsed time was %g seconds" % (end_time - start_time))
        start_time = time.time()
        dphi = p.d(x)
        end_time = time.time()
        print("Elapsed time was %g seconds" % (end_time - start_time))
        plt.plot(x, p(x))#, x, p.d(x))
        plt.show()

    def test_beta(self):
        return
        import scipy.stats
        a = 0.3
        b = 0.8
        rv = scipy.stats.beta(a, b)
        p = best.gpc.OrthogonalPolynomial(6, left=0, right=1, wf=rv.pdf)
        print 'At a single point: ', p(0.5)
        print p.d(0.5)
        print p([0.5, 0.3])
        print p.d([0.5, 0.3])
        quit()
        x = np.linspace(1e-4, 0.99, 10)
        print p(x)
        quit()
        #plt.plot(x, p(x))
        #plt.show()
        p = best.gpc.OrthogonalPolynomial(6, rv=rv)
        #plt.plot(x, p(x))
        #plt.show()
        rv = scipy.stats.expon()
        rv_cond = best.random.RandomVariableConditional(rv, (1, 2))
        print str(rv_cond)
        p = best.gpc.OrthogonalPolynomial(6, rv=rv_cond)
        x = np.linspace(1, 2, 100)
        plt.plot(x, p(x))
        plt.show()

    def test_product_basis(self):
        import time
        comp = (stats.beta(0.5, 0.5), stats.beta(0.5, 0.5))
        rv = best.random.RandomVectorIndependent(comp)
        print str(rv)
        prod = best.gpc.ProductBasis(degree=10, rv=rv)
        print str(prod)
        x = rv.rvs(size=10)
        print prod(x).shape
        x1 = np.linspace(1e-4, 0.99, 64)
        x2 = np.linspace(1e-4, 0.99, 64)
        X1, X2 = np.meshgrid(x1, x2)
        xx = np.vstack([X1.flatten(), X2.flatten()]).T
        z = rv.pdf(xx)
        Z = z.reshape((64, 64))
        #plt.contourf(X1, X2, np.log(Z).T)
        #plt.show()
        start_time = time.time()
        phi = prod(xx)
        end_time = time.time()
        print("Elapsed time was %g seconds" % (end_time - start_time))
        print phi.shape
        for j in range(phi.shape[1]):
            plt.contourf(X1, X2, phi[:, j].reshape((64, 64)))
            plt.show()


if __name__ == '__main__':
    unittest.main()
