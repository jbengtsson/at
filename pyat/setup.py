import glob
from io import open
import os
from os.path import abspath, basename, dirname, exists, join, splitext
import sys
import shutil
from setuptools import setup, Extension, find_packages

# Command line:
#   \rm -rf build
#   CC=g++-11 python3 setup.py <build|build_ext|install|--help> -i

# Numpy build dependency defined in pyproject.toml.
import numpy

print('\nPlatform:   ', sys.platform, '\nInstallation', end='')
if exists('/usr/local/include/omp.h'):
    print(' Homebrew ')
elif exists('/opt/local/include/libomp/omp.h'):
    print(' MacBook')


def select_omp():
    if exists('/usr/local/include/omp.h'):
        # Homebrew.
        return '-I/usr/local/include', '/usr/local/lib'
    elif exists('/opt/local/include/libomp/omp.h'):
        # MacBook.
        return '-I/opt/local/include/libomp', '/opt/local/lib/libomp'
    else:
        raise FileNotFoundError('\n'.join(('',
          'libomp.dylib must be installed with your favourite package manager:',
          '',
          'Use "$ brew install libomp"',
          'Or  "$ sudo port install libomp"',
          ''
        )))


here = abspath(dirname(__file__))
macros = [('PYAT', None)]
with_openMP = False


if True:
    # Clang.
    cflags = ['-std=gnu++14']
else:
    # GNU.
    cflags = [
        '-Wl,-no_compact_unwind', '-std=gnu++14', '-bundle',
        '-undefined dynamic_lookup', '-isysroot',
        '-I/usr/local/Cellar/gcc/11.2.0/include/c++/11.2.0/tr1',
        'limits.h']
#endif

if not sys.platform.startswith('win32'):
    cflags += ['-Wno-unused-function']

lflags = ['-L/usr/local/Cellar/armadillo/10.6.2/lib', '-larmadillo']


omp = os.environ.get('OPENMP', None)
if omp is None:
    omp_cflags = []
    omp_lflags = []
    omp_macros = []
else:
    # Get the location of an alternate OpenMP library
    # Example: OMP_MATLAB=$MATLABROOT/sys/os/glnx64
    omp_path = os.environ.get('OMP_MATLAB', None)
    # Get the threshold on the number of particles
    omp_threshold = int(os.environ.get('OMP_PARTICLE_THRESHOLD', 10))
    omp_macros = [('OMP_PARTICLE_THRESHOLD', omp_threshold)]
    if sys.platform.startswith('win'):
        omp_cflags = ['/openmp']
        omp_lflags = []
    elif sys.platform.startswith('darwin'):
        omp_inc, omp_lib = select_omp()
        omp_cflags = ['-Xpreprocessor', '-fopenmp', omp_inc]
        if omp_path is None:
            omp_lflags = ['-L' + omp_lib, '-lomp']
        else:
            omp_lflags = ['-L' + omp_path, '-Wl,-rpath,' + omp_path, '-liomp5']
    else:
        omp_cflags = ['-fopenmp']
        if omp_path is None:
            omp_lflags = ['-lgomp']
        else:
            omp_lflags = ['-L' + omp_path, '-Wl,-rpath,' + omp_path, '-liomp5']


# Get the long description from the README file
with open(join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


# It is easier to copy the integrator files into a directory inside pyat
# for packaging. However, we cannot always rely on this directory being
# inside a clone of at and having the integrator sources available, and
# this file is executed each time any setup.py command is run.
# It appears that only copying the files when they are available is
# sufficient.
at_source = abspath(join(here, 'at.cc'))
integrator_src_orig = abspath(join(here, '..', 'atintegrators'))
integrator_src = abspath(join(here, 'integrator-src'))
diffmatrix_source = abspath(
    join(here, '..', 'atmat', 'atphysics', 'Radiation')
)

if exists(integrator_src_orig):
    # Copy files into pyat for distribution.
#    source_files = glob.glob(join(integrator_src_orig, '*.[ch]'))
    source_files = glob.glob(join(integrator_src_orig, '*.h'))
    source_files.extend(glob.glob(join(integrator_src_orig, '*.cc')))
    source_files.extend(
        glob.glob(join(diffmatrix_source, 'diff_mat.cc'))
    )
    if not exists(integrator_src):
        os.makedirs(integrator_src)
    for f in source_files:
        shutil.copy2(f, integrator_src)

diffmatrix_method = join(integrator_src, 'diff_mat.cc')
elem_pass = join(integrator_src, 'ElemPass.cc')
at = Extension(
    'at.tracking.atpass',
    sources=[elem_pass, at_source, ],
    define_macros=macros + omp_macros,
    include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
    extra_compile_args=cflags + omp_cflags,
    extra_link_args=lflags + omp_lflags
)

diffmatrix = Extension(
    name='at.physics.diffmatrix',
    sources=[diffmatrix_method],
    include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
    define_macros=macros,
    extra_compile_args=cflags
)

setup(
    name='accelerator-toolbox',
    version='0.2.1',
    description='Accelerator Toolbox',
    long_description=long_description,
    author='The AT collaboration',
    author_email='atcollab-general@lists.sourceforge.net',
    url='https://github.com/atcollab/at',
    # Numpy 1.16.6 is the oldest version that builds with Python 3.9.
    install_requires=['numpy>=1.16.6', 'scipy>=0.16'],
    packages=find_packages(),
    ext_modules=[at, diffmatrix],
    zip_safe=False,
    python_requires='>=3.6.0',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ]
)
