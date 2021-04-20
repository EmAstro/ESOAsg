r"""Class to work on time cubes
"""

# import numpy as np
# from astropy.io import fits
# from astropy.io.fits.column import NUMPY2FITS
from astropy.table import Column, Table

# from ESOAsg import msgs
# from ESOAsg.core import fitsfiles
# from ESOAsg.ancillary import checks


from ESOAsg.datacontainers import datacontainer


class TimeCubes(datacontainer.DataContainer):
    r"""A class used to define and make simple operations on Time Cubes

    This allows to perform some basic tasks and to save files into  a format that is ESO Phase3 compliant.

    Attributes:

    Methods:

    """

    def __init__(self, primary_header=None, header=None, data=None, error=None, background=None, quality=None,
                 mask=None, others=None):
        r"""Instantiate the class LightCurves

        """
        self.primary_header = primary_header
        self.header = header
        self.data = data
        self.error = error
        self.background = background
        self.quality = quality
        self.mask = mask

        if others is not None:
            self.others = others
        else:
            self.others = Table()
        self._datatype = 'TimeCubes'
