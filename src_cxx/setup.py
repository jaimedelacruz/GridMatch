

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy
from distutils import sysconfig
#import numpy.distutils.intelccompiler
import numpy.distutils.ccompiler

#COMP = 'clang-6.0'
#COMPX = 'clang++-6.0'

COMP = 'gcc'
COMPX = 'g++'

os.environ["CC"] = COMP
os.environ["CXX"] = COMPX
#os.environ["F77"] = "gfortran"
#os.environ["FC"] = "gfortran"
from distutils import sysconfig
sysconfig.get_config_vars()['CFLAGS'] = ''
sysconfig.get_config_vars()['OPT'] = ''
sysconfig.get_config_vars()['PY_CFLAGS'] = ''
sysconfig.get_config_vars()['PY_CORE_CFLAGS'] = ''
sysconfig.get_config_vars()['CC'] = COMP #'clang'
sysconfig.get_config_vars()['CXX'] = COMPX #'clang++'
sysconfig.get_config_vars()['BASECFLAGS'] = ''
sysconfig.get_config_vars()['CCSHARED'] = ''
sysconfig.get_config_vars()['LDSHARED'] = COMP + " -shared"
sysconfig.get_config_vars()['CPP'] = COMPX
sysconfig.get_config_vars()['CPPFLAGS'] = ''
sysconfig.get_config_vars()['BLDSHARED'] = ''
sysconfig.get_config_vars()['CONFIGURE_LDFLAGS'] = ''
sysconfig.get_config_vars()['LDFLAGS'] = ''
sysconfig.get_config_vars()['PY_LDFLAGS'] = ''



comp_flags=['-O3','-std=c++17','-march=native','-fopenmp','-fPIC','-lstdc++']#,
root_dir = '/usr/'

setup(
    name = 'pyImtools',
    version = '1.0',
    author = 'Jaime de la Cruz Rodriguez (ISP-SU 2019)',
    # The ext modules interface the cpp code with the python one:
    ext_modules=[
        Extension("pyImtools",
            sources=["pyIMTOOLS.pyx", "libgrid.cpp","imtools.cc"], 
            include_dirs=["./", root_dir+"/include/",numpy.get_include()],
            language="c++",
            extra_compile_args=comp_flags,
                  extra_link_args=['-shared','-fopenmp'],#,'-lpython3.7m'],
                #      '-fPIC','-lgfortran','-fopenmp'],
            library_dirs=[root_dir+'/lib/','./'],#,'/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/'],#,'/Users/jaime/anaconda3/lib/'],
            libraries=[])
    ],
    cmdclass = {'build_ext': build_ext},
)
