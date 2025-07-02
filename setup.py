#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension
from wheel.bdist_wheel import bdist_wheel as native_bdist_wheel
import os
import subprocess
import platform
import re
import glob
try:
    from pip import main as pip_main
except Exception:
    from pip._internal import main as pip_main

try:
    from Cython.Build.Distutils import build_ext as native_build_ext
except ImportError:
    print('Suitable Cython unavailable, installing...')
    pip_main(['install', 'cython'])
    from Cython.Build.Distutils import build_ext as native_build_ext

try:
    import numpy as np
except ImportError:
    print('Suitable numpy unavailable, installing...')
    pip_main(['install', 'numpy'])
    import numpy as np


include_dirs = [np.get_include()]

# with open('README.rst') as readme_file:
#     readme = readme_file.read()

# with open('HISTORY.rst') as history_file:
#     history = history_file.read()

# requirements = [
#     'numpy>=1.16.0',
#     'pint>=0.7.0',
#     'xarray>=0.8.0',
#     'sympl==0.4.1',
#     'cython>=0.25',
#     'scipy>=0.18.1',
#     'climt==0.17.12',
# ]

# test_requirements = [
#     'pytest>=2.9.2',
#     'mock>=2.0.0',
# ]


# Create a custom build class to build libraries, and patch cython extensions
def build_libraries():

    # if os.environ.get('READTHEDOCS') == 'True':
    #     return

    curr_dir = os.getcwd()
    os.chdir(compiled_path)
    os.environ['PWD'] = compiled_path
    if subprocess.call(['make', 'CLIMT_ARCH='+operating_system]):
        raise RuntimeError('Library build failed, exiting')
    os.chdir(curr_dir)
    os.environ['PWD'] = curr_dir


# Custom build class
class gfs_build_ext(native_build_ext):

    def run(self):
        build_libraries()
        native_build_ext.run(self)


# Custom bdist_wheel class
class gfs_bdist_wheel(native_bdist_wheel):

    def run(self):
        self.run_command('build')
        native_bdist_wheel.run(self)

setup(name='pyTRACK',
      version='1.0',
      description='python wrapper for TRACK',
      author='Ai33L',
      packages=['pyTRACK'],
     )
