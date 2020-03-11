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

setup(
    name='mamonca',
    version='0.0.1',
    description='mamonca - an integrated development environment (IDE) for computational materials science.',
    url='https://github.com/samwaseda/sams-studio/tree/master/MC',
    author='Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department',
    author_email='waseda@mpie.de',
    license='BSD',
    cmdclass = {'build_ext': my_build_ext}, ext_modules=cythonize("mc.pyx", language_level="3")
)
#setup(ext_modules=cythonize("mc.pyx", language_level="3"))
