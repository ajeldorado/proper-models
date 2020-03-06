
import os
import sys
import platform
import ez_setup
ez_setup.use_setuptools()

from setuptools import find_packages, setup, Extension

copy_args = sys.argv[1:]
copy_args.append('--user')
 
ext_modules = []

setup(
      name="wfirst_phaseb_proper",
      version = "1.7",
      packages=find_packages(),

      install_requires = ['numpy>=1.8', 'scipy>=0.19', 'astropy>=1.3', 'PyPROPER>=3.1.5;python_version<"3.0"', 
        'PyPROPER3>=3.1.5;python_version>="3.0"'],

      package_data = {
        '': ['*.*']
      },

      script_args = copy_args,

      zip_safe = False, 

      # Metadata for upload to PyPI
      author="John Krist",
      author_email = "john.krist@jpl.nasa.gov",
      description="WFIRST coronagraph Phase B PROPER prescription",
      license = "BSD",
      platforms=["any"],
      url="",
      ext_modules = ext_modules
)
