=================
Installing ESOAsg
=================

This document describes how to install ESOAsg.

Installing Dependencies
=======================

There are a few packages that need to be installed before running the package.

We highly recommend that you use Anaconda for the majority of these installations.

Detailed installation instructions are presented below:

python and dependencies
-----------------------

* `python <http://www.python.org/>`_ version 3.7 or later
* `astropy <https://www.astropy.org/>`_ -- version 4.0 or later
* `astroquery <https://astroquery.readthedocs.io/en/latest/>`_ -- version 0.4 or later
* `IPython <https://ipython.org>`_ -- version 7.13 or later
* `matplotlib <https://matplotlib.org/>`_ -- version 3.2 or later
* `numpy <http://www.numpy.org/>`_ -- version 1.18.0 or later
* `pyvo <https://pypi.org/project/pyvo/>`_ -- version 1.0 or later

https://packaging.python.org/

packaging>=19.0

photutils>=0.7.2
pytest>=3.0.7
>=1.0


git clone
---------

To install the package via GitHub run::

    git clone https://github.com/EmAstro/ESOAsg.git

And we then recommend you install with::

    python setup.py develop
