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


def from_element_to_list(element, element_type=str):
    r"""Given an element it returns a list containing the element

    It also checks all the elements in the list have the same type defined by `element_type`

    Args:
        element (any): element that will be put in the string
        element_type (any): type of the element that should be contained in the list

    Returns:
        list: list containing `element`

    """
    if element is None:
        return None
    elif isinstance(element, list):
        for element_in_list in element:
            assert isinstance(element_in_list, element_type), r'{} must be a {}'.format(element_in_list, element_type)
        return element
    elif isinstance(element, element_type):
        return [element]
    else:
        msgs.error('Not valid type for: {}'.format(element))
    return


def from_number_to_string(number):
    r"""Given an int or a float it returns a string

    If the input is a int, it is first converted to float and then to string

    Args:
        number (any): `int` or `float` that needs to be transformed into a string

    Returns:
        str: same of number but as a `str`

    """
    if number is None:
        return None
    elif isinstance(number, str):
        return number
    elif isinstance(number, int):
        return str(float(number))
    elif isinstance(number, float):
        return str(number)
    else:
        msgs.error('The value entered is not a string or a number. Type: {}'.format(type(number)))
        return


def from_bytes_to_string(input_in_bytes):
    r"""Given an input in `bytes` return it the corresponding `str`

    This is mainly to deal with the fact that TAP queries return a list in bytes format that might be annoying. If
    the input is not `bytes` nothing is changed.

    Args:
        input_in_bytes (bytes): input in bytes

    Returns:
        str: output converted into a string

    """
    if isinstance(input_in_bytes, bytes):
        output_as_str = np.str(input_in_bytes.decode("utf-8"))
    elif isinstance(input_in_bytes, np.ndarray) and len(input_in_bytes.shape) == 1:
        output_as_str = np.empty_like(input_in_bytes)
        for idx in np.arange(0, np.size(output_as_str)):
            if isinstance(input_in_bytes[idx], bytes):
                output_as_str[idx] = np.str(input_in_bytes[idx].decode("utf-8"))
            else:
                output_as_str[idx] = input_in_bytes[idx]
    else:
        output_as_str = input_in_bytes
    return output_as_str


def remove_non_ascii(text_string):
    r"""Replace non ascii characters from a string

    Args:
        text_string (`str`):
            input string from which the non ascii characters will be removed

    Returns:
        text_string_cleaned ('str'):
            string from which the non ASCII characters have been removed.

    """
    text_string_cleaned = "".join(character for character in text_string if 31 < ord(character) < 123)
    text_string_cleaned.replace('\\', '').strip()
    return text_string_cleaned


def check_disk_space(min_disk_space=float(default.get_value('min_disk_space'))):
    r"""
    Given a limit in GB in the variable min_disk_space the macro returns `True` if there is enough space and rises
    an error otherwise.

    Args:
        min_disk_space (`float`):
            Size of free space on disk required

    Returns:
        enough_space (`bool`):
            True if there is enough space on disk
    """
    total, used, free = shutil.disk_usage("./")
    total = total / (1024 ** 3)
    used = used / (1024 ** 3)
    free = free / (1024 ** 3)
    msgs.info('Your disk has: Total: {0:.2f} GB, Used: {0:.2f} GB, Free: {0:.2f} GB'.format(total, used, free))
    if free > min_disk_space:
        enough_space = np.bool(1)
    else:
        enough_space = np.bool(0)
        msgs.error('Not enough space on disk')
    return enough_space


def connection_to_website(url, timeout=1):  # written by Ema 05.03.2020
    r"""Check there is an active connection to a website

    Args:
        url (`str`):
            link to the website you want to check
        timeout (`int`, `float`):
            timeout waiting for the website to respond

    Returns:
         is_active (`boolean`):
            `True` if there is an active connection, `False` and error raised if not.

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


def fits_file_is_valid(fits_file, verify_fits=False, overwrite=False):  # Written by Ema 05.03.2020
    r"""Check if a file exists and has a valid extension

    The option `verify_fits` checks the header of the fits file using `astropy.io.fits.verify`

    Args:
        fits_file (`str`):
            fits file you would like to check
        verify_fits (`bool`):
            if set to `True`, it will verify that the fits file is complaint to the FITS standard.
        overwrite (`bool`):
            if `True`, overwrite the input fits file with the header corrections from `verify_fits`

    Returns:
        is_fits (`boolean`):
            `True` if exists `False` and warning raised if not.

    """
    is_fits = True
    # Checks if it is a string
    assert isinstance(fits_file, str), 'input fits needs to be a string'
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


def check_checksums(hdul):
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


def image2d_is_valid(image2d):  # Written by Ema 12.03.2020
    r"""Check if a 2D image is valid

    Args:
        image2d (np.array):
            image that you would like to check

    Returns:
        is_image2d (`boolean`):
            `True` if a valid 2D image `False` and error raised if not.
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


'''
def single_value_to_list(single_value):
    """This is useful to transform single strings, integers, floats into python lists. Can be used when the
    input to a function is given by the `parse_arguments()`

    Args:
        single_value: (str, float, int):
            This is the argument you want to have as a list

    Returns:
        list_value: (list)
            List containing the values given in input
    """
    print(single_value)
    print(type(single_value))
    if type(single_value) is 'list':
        list_value = single_value
    elif isinstance(single_value, str):
        list_value = [single_value]
    else:
        print('Noting')
        list_value = single_value
    print(' ')
    print(list_value)
    print(type(list_value))
    return list_value



# ToDo:
def check_instrument(instrument):
    """Given an instrument name, it checks if it is
    a valid entry

    Args:
        instrument (str):
            Instrument name you want to check

    Returns:
        is_instrument (np:bool):
            True if it is a valid instrument
    """
    instrument_list = ['MUSE']
    if instrument in instrument_list:
        is_instrument = True
    else:
        is_instrument = False
        msgs.error('Wrong instrument name')
    return is_instrument
'''
