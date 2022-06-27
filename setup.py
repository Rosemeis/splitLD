from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [Extension(
				"splitLD.shared_cy",
				["splitLD/shared_cy.pyx"],
				extra_compile_args=['-fopenmp', '-g0'],
				extra_link_args=['-fopenmp'],
				include_dirs=[numpy.get_include()]
			)]

setup(
	name="splitLD",
	version="1.0",
	description="Optimal linkage disequilibrium splitting",
	author="Jonas Meisner",
	packages=["splitLD"],
	entry_points={
		"console_scripts": ["splitLD=splitLD.main:main"]
	},
	python_requires=">=3.6",
	ext_modules=cythonize(extensions),
	include_dirs=[numpy.get_include()]
)
