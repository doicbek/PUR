from numpy.distutils.core import setup, Extension
import setuptools
from Cython.Build import cythonize
from numpy import get_include
import os

nmpy_inc = get_include()

setup(name='pur',
	version="0.0.1",
	description='lightweight pure pseudo-Cl estimator for HEALPix',
	author='Ari Cukierman, Dominic Beck',
	author_email='dobeck@stanford.edu',
	packages=['pur'],
	ext_modules=cythonize([Extension("pur.mcm", sources=["pur/mcm.pyx"], include_dirs=[nmpy_inc], extra_compile_args=["-O3"], language="c++")]),
	install_requires=['numpy','healpy'],
	license='GPLv2'
)