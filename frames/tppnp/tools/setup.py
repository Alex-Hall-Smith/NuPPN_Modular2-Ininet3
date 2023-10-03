from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(name = 'smooth',
      ext_modules = cythonize([Extension("smooth",
                     sources=["smooth.pyx"],
                     libraries=["m"],
                     extra_compile_args=["-O3","-std=c99"]
                     )
                 ]
           )
      )
