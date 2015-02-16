#!/usr/bin/env python
try:
    from setuptools import setup
except:
    from distutils.core import setup

# -------------------------------------------------------------------- #
# Run Setup
setup(name='rasmlib',
      version='0.0.0',
      description='Tools library for the Regional Arctic System Model',
      long_description='',
      author='Joe Hamman',
      author_email='jhamman@hydro.washington.edu',
      install_requires=['scipy >= 0.13', 'numpy >= 1.8',
                        'netCDF4 >= 1.0.6', 'matplotlib >= 1.3.1',
                        'nco >= 0.0.2'],
      url='https://github.com/jhamman/rasmlib',
      packages=['rasmlib', 'rasmlib.analysis', 'rasmlib.post_processing'],
      scripts=['scripts/rasm_analysis', 'scripts/rasm_post_process'])
# -------------------------------------------------------------------- #
