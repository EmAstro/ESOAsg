"""
Class to work on images
"""

from astropy.io import fits

from ESOAsg import msgs
from ESOAsg.ancillary import checks


class Image:
    r"""A class used to define and make simple operations on images

    This allows to perform some basic checks on a 2D image.

    Attributes:

    Methods:

    """

    def __init__(self, data=None, errors=None, background=None, quality=None, bad_pixels_mask=None, sources_mask=None):
        r"""Instantiate the class image

       """

    def from_fits(self, fits_file=None):
        r"""

        Args:
            fits_file (`str`):
                fits file name where the image is stored

        Returns:
            image (`Image')
                a filled Image
        """
        checks
        hdul = fits.open(fits_file)
        
