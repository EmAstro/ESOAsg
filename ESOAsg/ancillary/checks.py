r"""Module that performs some useful and basic checks and transformations
"""

# import sys
import numpy as np
import shutil
import urllib
import os.path

from astropy.io import fits

# ESOAsg imports
from ESOAsg import msgs
from ESOAsg import default


def check_disk_space(min_disk_space=float(default.get_value('min_disk_space'))) -> bool:
    r"""Check that there is enough space on the location where the code is running

    Given a disk space limit in GB, the macro returns `True` if the disk where the code is running has more free GB
    than the given limit.

    .. warning::
        The current implementation checks the disk where the code is running (i.e., from the directory: `./`).
        This may cause some troubles with shared disks.

    Args:
        min_disk_space (float): Size of free space on disk required

    Returns:
        bool: `True` if there is enough space on disk

    """
    total, used, free = shutil.disk_usage("./")
    total = total / (1024. ** 3.)
    used = used / (1024. ** 3.)
    free = free / (1024. ** 3.)
    msgs.info('Your disk has:')
    msgs.info('Total: {0:.2f} GB, Used: {0:.2f} GB, Free: {0:.2f} GB'.format(total, used, free))
    if free > min_disk_space:
        enough_space = True
    else:
        enough_space = False
        msgs.warning('Not enough space on disk')
    return enough_space


def connection_to_website(url, timeout=1.) -> bool:  # written by Ema 05.03.2020
    r"""Check there is an active connection to a website

    Args:
        url (str): link to the website you want to check
        timeout (float): timeout waiting for the website to respond

    Returns:
         bool: `True` if there is an active connection, `False` and error raised if not.

    """
    # Checks for url
    assert isinstance(url, str), 'The url needs to be a string'
    if url.startswith('www'):
        url_clean = 'http://' + url
        msgs.warning('Modifying url to: {}'.format(url_clean))
    else:
        url_clean = url

    request = urllib.request.Request(url_clean)
    try:
        urllib.request.urlopen(request, timeout=timeout)
    except urllib.error.HTTPError as err:
        msgs.warning('HTTP Error: {}'.format(err.code))
        is_active = False
    except urllib.error.URLError as err:
        msgs.warning('URL Error: {}'.format(err.reason))
        is_active = False
    else:
        is_active = True
    return is_active


def fits_file_is_valid(fits_file, verify_fits=False, overwrite=False) -> bool:  # Written by Ema 05.03.2020
    r"""Check if a file exists and has a valid extension

    The option `verify_fits` checks the header of the fits file using `astropy.io.fits.verify`

    Args:
        fits_file (str): fits file you would like to check
        verify_fits (bool): if set to `True`, it will verify that the fits file is complaint to the FITS standard.
        overwrite (bool): if `True`, overwrite the input fits file with the header corrections from `verify_fits`

    Returns:
        bool: `True` if exists `False` and warning raised if not.

    """
    is_fits = True
    # Checks if it is a string
    assert isinstance(fits_file, str), 'input `fits_file` needs to be a string'
    # Check for ending
    if not fits_file.endswith('.fits') and not fits_file.endswith('.fits.fz') and not fits_file.endswith('.fits.gz'):
        msgs.warning('File: {} does not end with `fits` or `fits.fz` or `fits.gz`'.format(fits_file))
        is_fits = False
    # Check for existence
    if not os.path.exists(fits_file):
        msgs.warning('File: {} does not exists'.format(fits_file))
        is_fits = False
    # Check for compliance with FITS standard
    if verify_fits:
        if overwrite:
            hdul = fits.open(fits_file, mode='update', checksum=False)
            if not check_checksums(hdul):
                is_fits = False
            hdul.flush(output_verify='fix+warn', verbose=True)
            hdul.writeto(fits_file, checksum=True, overwrite=True)
            msgs.info('File checked and rewritten')
        else:
            hdul = fits.open(fits_file, mode='readonly', checksum=True)
            if not check_checksums(hdul):
                is_fits = False
            hdul.verify('fix+warn')
        hdul.close()
    else:
        if overwrite:
            msgs.error('The option overwrite works only if verify_fits = True')
    return is_fits


def check_checksums(hdul) -> bool:
    r"""Test if the `datasum` and `checksum` keywords in a `HDUList` are present and up-to-date

    Args:
        hdul (:py:obj:`astropy.io.fits.hdu.hdulist.HDUList`): list of `astropy` HDUs to be checked

    Returns:
        bool: `True` all the HDUs in the input HDUL have the correct `datasum` and `checksum`

    """
    is_good_checksum = True
    for hdu in hdul:
        checks_for_checksum = hdu.verify_checksum()
        checks_for_datasum = hdu.verify_datasum()
        if checks_for_checksum == 0:
            msgs.warning('Checksum not valid')
            is_good_checksum = False
        if checks_for_checksum == 2:
            msgs.warning('Checksum not present')
            is_good_checksum = False
        if checks_for_datasum == 0:
            msgs.warning('Datasum not valid')
            is_good_checksum = False
        if checks_for_datasum == 2:
            msgs.warning('Datasum not present')
            is_good_checksum = False
    return is_good_checksum


def image2d_is_valid(image2d) -> bool:  # Written by Ema 12.03.2020
    r"""Check if a 2D image is valid

    Args:
        image2d (obj:`numpy.ndarray`): image that you would like to check

    Returns:
        bool: `True` if a valid 2D image `False` and error raised if not.
    """
    is_image2d = True

    # Checks if it is a numpy array
    assert isinstance(image2d, np.ndarray), 'The image is not a `numpy array`'
    # Check for dimensions
    if not image2d.ndim == 2:
        msgs.warning('The image is not two dimensional (N. of dimension={})'.format(image2d.ndim))
        is_image2d = False

    return is_image2d


def table_is_valid(table):  # Written by Ema 08.04.2020
    r"""Check if a table is valid

    Args:
        table (`fits.BinTableHDU` or `fits.TableHDU`):
            table that you would like to check

    Returns:
        is_table (`boolean`):
            `True` if a valid table format `False` and error raised if not.
    """
    is_table = True

    # Checks if it is an astropy table
    assert isinstance(table,
                      (fits.BinTableHDU, fits.TableHDU)), 'The table is not a `fits.BinTableHDU` or a `fits.TableHDU`'

    return is_table


def header_is_valid(header):
    r"""Check if an header is valid

    """
    is_header = True

    # Check if is a fits.header
    assert isinstance(header, fits.Header), 'The header is not an instance of `astropy.fits.io.header`'
    if len(header) == 0:
        msgs.warning('Empty Header')
        is_header = False

    return is_header
