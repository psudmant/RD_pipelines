#cython: language_level=3
#cython: linetrace=False

import os
import numpy as np

from distutils.core import setup, Extension 
from Cython.Distutils import build_ext

numpy_path = os.path.basename(np.__file__) + "/core/include/"


ext_modules = [
               Extension("calc_gc_depth",
                         ["calc_gc_depth.pyx"],
                         language="c",
                         include_dirs = [np.get_include()]
                        )
              ]

setup( 
	name = 'calc gc depth app',
	cmdclass = {'build_ext': build_ext},
	ext_modules = ext_modules,
	)
