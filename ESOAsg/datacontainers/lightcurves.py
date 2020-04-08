r"""Class to work on light curves
"""


import numpy as np
from astropy.io import fits


from ESOAsg import msgs
from ESOAsg.core import fitsfiles
from ESOAsg.ancillary import checks


class LightCurves:
    r"""A class used to define and make simple operations on time series

    This allows to perform some basic tasks on a time series and to save it in a format that is ESO Phase3
    compliant.

    Attributes:

    Methods:

    """

    def __init__(self, header=None, time=None, time_bin=None, flux=None, mag=None, error=None,
                 background=None, quality=None, others=None):
        r"""Instantiate the class LightCurves


        Each field of the BINTABLE shall be further described in the extension header. Mandatory fields shall be:
        * time
        * flux (or mag)
        * error
        in that particular order.
        Additional fields may be added.



        """
        self.header = header
        self.time = time
        self.time_bin = time_bin
        self.flux = flux
        self.mag = mag
        self.error = error
        self.background = background
        self.quality = quality
        self.others = others

    def load_from_table(self, table, where_time='TIME', where_time_bin='TIME_BIN', where_flux='FLUX', where_mag='MAG',
                        where_error='ERROR', where_background='BACKGROUND', where_quality='QUAL'):
        r"""Given an astropy table put it in a LightCurves object
        """
        for attribute, value in zip(self.__dict__.keys(), self.__dict__.values()):
            if attribute is not 'header':
                where_attribute = vars()['where_' + attribute]
                print(where_attribute)
                from IPython import embed
                embed()
                if where_attribute in table.colnames:
                    setattr(self, attribute, table.data[where_attribute])

        return True

