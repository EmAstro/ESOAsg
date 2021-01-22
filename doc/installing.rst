=================
Installing ESOAsg
=================

This section describes how to install ESOAsg.

Installing Dependencies
=======================

There are a few packages that need to be installed before running the package.

We highly recommend that you use Anaconda for the majority of these installations.

Detailed installation instructions are presented below:

python and dependencies
-----------------------

ESOAsg runs with `python <http://www.python.org/>`_ 3.7 and with the following dependencies:

* `python <http://www.python.org/>`_ -- version 3.7 or later
* `astropy <https://www.astropy.org/>`_ -- version 4.2 or later
* `astroquery <https://astroquery.readthedocs.io/en/latest/>`_ -- version 0.4 or later
* `IPython <https://ipython.org>`_ -- version 7.13 or later
* `matplotlib <https://matplotlib.org/>`_ -- version 3.3 or later
* `numpy <http://www.numpy.org/>`_ -- version 1.19 or later
* `packaging <https://packaging.python.org/>`_ -- version 20.4 or later
* `pytest <https://docs.pytest.org/>`_ version 6.2 or later
* `pyvo <https://pypi.org/project/pyvo/>`_ -- version 1.1 or later
* `requests <https://requests.readthedocs.io/>`_ -- version 2.25 or later

The following packages are required to work on Ligo gravitational wave events, see for instance `this notebook <https://github.com/EmAstro/ESOAsg/blob/master/doc/notebooks/HOWTO_getDataFromGWContours.ipynb>`_:

* `healpy <https://healpy.readthedocs.io/>`_ -- version 1.14 or later
* `ligo.skymap <https://lscsoft.docs.ligo.org/ligo.skymap/>`_ -- version 0.5 or later

If you are using Anaconda, you can check the presence of these packages with::

    conda list "^python$|astropy$|astroquery$|IPython|matplotlib|numpy|packaging|pytest|pyvo|requests|healpy|ligo.skymap"

If the packages have been installed, this command should print out all the packages and their version numbers.

If any of the packages are out of date, they can be updated with a command like::

    conda update astropy

You might have troubles in installing `pyvo <https://pypi.org/project/pyvo/>`_.
In case, you can get it from GitHub running in the directory where you want to install it::

    git clone https://github.com/astropy/pyvo.git
    cd pyvo
    python setup.py install

git clone
---------

To install the package via GitHub run::

    git clone https://github.com/EmAstro/ESOAsg.git

And, given that the packages is still work in progress and you may want to updated on-the-fly, we then recommend to install it with the `develop` option::

    cd ESOAsg
    python setup.py develop

Testing the Installation
========================

In order to assess whether ESOAsg has been properly installed, we suggest you run the following tests:

1. Run the default tests
------------------------

In the directory where ESOAsg is installed run::

    pytest

You should see that all the current test passed.

2. Ensure that the scripts work
-------------------------------

Go to a directory outside of the ESOAsg directory, then type mod_header::

    cd
    mod_header -h

