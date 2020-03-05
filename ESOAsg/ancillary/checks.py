"""
Module that performs some useful and basic checks
"""

# import sys
import numpy as np
import shutil
import urllib
import os.path

# ESOAsg imports
from ESOAsg import msgs
from ESOAsg import default


def check_disk_space(min_disk_space=np.float32(default.get_value('min_disk_space'))):
    r"""
    Given a limit in GB in the variable min_disk_space the macro returns `True` if there is enough space and rises
    an error otherwise.

    Args:
        min_disk_space (`numpy.float32`):
            Size of free space on disk required

    Returns:
        enough_space (`numpy.bool`):
            True if there is enough space on disk
    """
    total, used, free = shutil.disk_usage("./")
    total = total / (1024**3)
    used = used / (1024**3)
    free = free / (1024**3)
    msgs.info('Your disk has: Total: {0:.2f} GB, Used: {0:.2f} GB, Free: {0:.2f} GB'.format(total, used, free))
    if free > min_disk_space:
        enough_space = np.bool(1)
    else:
        enough_space = np.bool(0)
        msgs.error('Not enough space on disk')
    return enough_space


def connection_to_website(url, timeout=1):   # written by Ema 05.03.2020
    r"""Check there is an active connection to a website

    Args:
        url (`str`):
            link to the website you want to check
        timeout (`int`, `float`):
            timeout waiting for the website to respond

    Returns:
        `boolean`:
            `True` if there is an active connection, `False` and error raised if not.

    """
    # Checks for url
    assert isinstance(url, str), 'The url needs to be a string'
    if url.startswith('www'):
        url_clean = 'http://'+url
        msgs.warning('Modifying url to: {}'.format(url_clean))
    else:
        url_clean = url

    request = urllib.request.Request(url_clean)
    try:
        urllib.request.urlopen(request, timeout=timeout)
    except urllib.error.HTTPError as err:
        msgs.warning('HTTP Error: {}'.format(err.code))
        return False
    except urllib.error.URLError as err:
        msgs.warning('URL Error: {}'.format(err.reason))
        return False
    else:
        return True


def fits_file_is_valid(fits_file):  # Written by Ema 05.03.2020
    r"""Check if a file exists and has a valid extension

    Args:
        fits_file (`str`):
            fits file you would like to check

    Returns:
        boolean`:
            `True` if exists `False` and error raised if not.

    """

    # Checks for url
    assert isinstance(fits_file, str), 'input fits needs to be a string'
    # Check for ending
    if not fits_file.endswith('.fits') and not fits_file.endswith('.fits.fz'):
        msgs.warning('File: {} does not end with `.fits` or .`fits.fz`'.format(fits_file))
        return False
    # Check for existence
    if os.path.exists(fits_file):
        return True
    else:
        msgs.warning('File: {} does not exists'.format(fits_file))
        return False


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
    instrument_list = ['MUSE', 'SINFONI']
    if instrument in instrument_list:
        is_instrument = True
    else:
        is_instrument = False
        msgs.error('Wrong instrument name')
    return is_instrument
'''
