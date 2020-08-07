from distutils.core import setup
from Cython.Build import cythonize

from distutils.command.build_ext import build_ext
from distutils.sysconfig import customize_compiler

class my_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        build_ext.build_extensions(self)

setup(cmdclass = {'build_ext': my_build_ext},
      ext_modules=cythonize("tree.pyx", language_level="3"),
      include_dirs=['/u/samsstud/local/include/EIGEN'])
