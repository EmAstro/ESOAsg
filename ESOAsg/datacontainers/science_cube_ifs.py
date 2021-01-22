r"""Class to work on SCIENCE.CUBE.IFS
"""

# import numpy as np
# from astropy.io import fits
# from astropy.io.fits.column import NUMPY2FITS
from astropy.table import Column, Table

# from ESOAsg import msgs
# from ESOAsg.core import fitsfiles
# from ESOAsg.ancillary import checks
from ESOAsg.datacontainers import datacontainer


class ScienceCubeIfs(datacontainer.DataContainer):
    r"""A class used to define and make simple operations on data of type PRODCATG=`SCIENCE.CUBE.IFS`

    This is a child of :obj:`ESOAsg.datacontainers.DataContainer`.

    Please refer to the `ESO Science Data Products Standard <https://www.eso.org/sci/observing/phase3/p3sdpstd.pdf>`_
    for a full description of the format.

    Attributes:

    Methods:

    """

    _prodcatg_type = 'SCIENCE.CUBE.IFS'

    def __init__(self, primary_header=None, scidata_hdu=None, errdata_hdu=None, qualdata_hdu=None,
                 others_hdu=None):
        r"""Instantiate the class ScienceCubeIfs
        """
        super().__init__(prodcatg_type=self._prodcatg_type)
        self.primary_header = primary_header
        self.scidata = scidata_hdu
        self.errdata = errdata_hdu
        self.qualdata = qualdata_hdu
        self.others = others_hdu


