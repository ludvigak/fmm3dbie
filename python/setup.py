import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from sys import platform

pkg_name = "solvers3dpy"

list_files=[]
list_files.append('../src/helm_wrappers/helm_comb_dir.f')
list_files.append('../src/common/sparse_reps.f')
list_files.append('../src/surface_routs/surf_routs.f90')

FAST_KER = os.getenv('FAST_KER')
FFLAGS = os.getenv('FFLAGS')

FLIBS=[]
FLIBS.append('-lm')
FLIBS.append('-lgomp')
FLIBS.append('-lfmm3d')
FLIBS.append('-L../lib')
if platform == "linux" or platform == "linux2":
    FLIBS.append('-lopenblas')
FLIBS.append('../lib-static/libsolvers3d.a')

FFLAGS = FFLAGS.rstrip().split(' ')
FFLAGS = list(filter(None,FFLAGS))

helm = []
com = []
surf = []
helm.append('helm_comb_dir_fds_mem')
helm.append('helm_comb_dir_fds_init')
helm.append('helm_comb_dir_fds_matgen')
helm.append('lpcomp_helm_comb_dir')
com.append('conv_to_csc')
surf.append('surf_vals_to_coefs')
surf.append('get_qwts')
surf.append('get_patch_id_uvs')

ext_helm = Extension(
    name='helm3d_dir',
    sources=list_files,
    f2py_options=['only:']+helm+com+surf+[':'],
    extra_f77_compile_args=FFLAGS,
    extra_f90_compile_args=FFLAGS,
    extra_link_args=FLIBS
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    version="0.1.0",
    author="Manas Rachh",
    author_email="mrachh@flatironinstitute.org",
    description="This pacakge contains basic routines Helmholtz dirichlet fast direct solver",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
        "pytest"
    ],
    ext_modules=[ext_helm],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )    
)
