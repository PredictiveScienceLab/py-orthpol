#!/usr/bin/env python


from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import os
import glob


setup(name='Orthogonal Polynomials in Python',
      author='Ilias Bilionis',
      version='0.0',
      ext_modules=[Extension('orthpol._orthpol',
                             glob.glob(os.path.join('src',
                                                    '*.f')),
                             libraries=['python2.7', 'pthread', 'tcl8.5'],
                             extra_link_args=['-shared'])],
      packages=['orthpol'])
