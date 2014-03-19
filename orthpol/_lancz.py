"""

Author:
    Ilias Bilionis

Date:
    7/25/2013
"""


__all__ = ['lancz']


import numpy as np
import _orthpol as orthpol


def lancz(x, w, n):
    """The Lanczos procedure for constructing the recurcive formula
    for orthogonal polynomials.

    Wrapper from ORTHPOL.
    """
    if x.dtype == 'float32':
        func = orthpol.slancz
    else:
        func = orthpol.dlancz
    alpha, beta, ierr = func(n, x, w)
    assert ierr == 0
    return alpha, beta


if __name__ == '__main__':
    x = np.linspace(0., 1, 100)
    w = np.ones(100)
    print x, w
    alpha, beta = lancz(x, w, 10)
    print alpha
    print beta
