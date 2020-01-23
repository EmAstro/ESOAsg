"""
Create header class
"""

from astropy.io import fits
from os import path

from ESOAsg import msgs


class Header:
    r"""
    Instantiate an header.
    This makes large use of `astropy.fits.PrimaryHDU()`
    """

    def __init__(self):
        self.empty = fits.PrimaryHDU().header

    def from_file(self, fits_name=None):
        r"""
        Load an header with the information from a fits file

        Args:
            fits_name (`str`):
                fits file name
        """

        if fits_name is None:
            msgs.warning('No file selected, returning empty header')
            self.empty
        elif not path.exists(fits_name):
            msgs.warning('File not exist, returning empty header')
            self.empty
        else:
            hdul = fits.open(fits_name)
            print(hdul.header)
