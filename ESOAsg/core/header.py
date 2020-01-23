"""
Module to hack the header of a file
"""

from astropy.io import fits
from os import path

from ESOAsg import msgs


def from_file(fits_name=None, which_hdu=0):
    r"""
    Load an header with the information from a fits file

    Args:
        fits_name (`str`):
            fits file name
        which_hdu (`numpy.int`):
            select from which HDU you are getting the header

    Returns:
         header (`fits.header`):
             the header corresponding to `which_hdu` from `fits_name`
    """
    if fits_name is None:
        msgs.error('No file selected.')
    elif not path.exists(fits_name):
        msgs.warning('File does not exist.')
    else:
        hdul = fits.open(fits_name)
        msgs.info('The fits file contains {} HDUs'.format(len(hdul)))
        return hdul[which_hdu].header



