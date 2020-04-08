from distutils.core import setup
from distutils.extension import Extension

import Cython.Compiler.Options
from Cython.Build import cythonize

Cython.Compiler.Options.annotate = True

ext_modules = [
    Extension("genometargeting._compute",
              sources=["genometargeting/_compute.pyx"],
              libraries=["m"],  # Unix-like specific
              extra_compile_args=["-O3", "-fopenmp"],
              extra_link_args=['-fopenmp'],
              )
]

setup(name="Compute",
      ext_modules=cythonize(ext_modules),
      )
