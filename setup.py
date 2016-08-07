import distutils
from distutils.core import setup, Extension, Command
import os
import numpy

data_files=[]

ext=Extension("cosmology._cosmolib",
              ["cosmology/cosmolib_pywrap.c", "cosmology/cosmolib.c"])

include_dirs=[numpy.get_include()]

setup(name="cosmology",
      packages=['cosmology'],
      version="1.1.0",
      data_files=data_files,
      ext_modules=[ext],
      include_dirs=include_dirs)




