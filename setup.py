# !usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys
import os
import glob

from setuptools import setup, find_packages

def get_scripts():
    """ Grab all the scripts in the bin directory.  """
    scripts = []
    if os.path.isdir('bin'):
        scripts = [ fname for fname in glob.glob(os.path.join('bin', '*'))
                                if not os.path.basename(fname).endswith('.rst') ]
    return scripts


def get_requirements():
    """ Get the requirements from a system file.  """
    name = 'ESOAsg/requirements.txt'

    requirements_file = os.path.join(os.path.dirname(__file__), name)
    install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                        if not line.strip().startswith('#') and line.strip() != '']
    return install_requires


NAME = 'ESOAsg'
VERSION = '0.00'
AUTHOR = 'Ema'


def run_setup(data_files, scripts, packages, install_requires):

    # TODO: Are any/all of the *'d keyword arguments needed? I.e., what
    # are the default values?

    setup(name=NAME,
          provides=NAME,                                                # *
          version=VERSION,
          license='TBD',
          description='ESOAsg Usefull Tools',
          long_description=open('README.md').read(),
          author='PypeIt Collaboration',
          author_email='pypeit@ucolick.org',
          keywords='pypeit PypeIt astronomy Keck UCO Lick data reduction',
          url='https://github.com/pypeit/PypeIt',
          packages=packages,
          package_data={'pypeit': data_files, '': ['*.rst', '*.txt']},
          include_package_data=True,
          scripts=scripts,
          install_requires=install_requires,
          requires=[ 'Python (>3.6.0)' ],                               # *
          zip_safe=False,                                               # *
          use_2to3=False,                                               # *
          setup_requires=[ 'pytest-runner' ],
          tests_require=[ 'pytest' ],
          classifiers=[
              'Development Status :: 2 - Pre-Alpha',
              'Intended Audience :: Science/Research',
              'License :: TBD',
              'Natural Language :: English',
              'Operating System :: OS Independent',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3.6',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Topic :: Software Development :: Libraries :: Python Modules',
              'Topic :: Software Development :: User Interfaces'
          ],
          )

#-----------------------------------------------------------------------
if __name__ == '__main__':

    # Compile the data files to include
    data_files = get_data_files()
    # Compile the scripts in the bin/ directory
    scripts = get_scripts()
    # Get the packages to include
    packages = find_packages()
    # Collate the dependencies based on the system text file
    install_requires = get_requirements()
    install_requires = []  # Remove this line to enforce actual installation
    # Run setup from setuptools
run_setup(data_files, scripts, packages, install_requires)
