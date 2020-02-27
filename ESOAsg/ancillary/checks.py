"""
Module that performs some useful and basic checks
"""

# import sys
import shutil
import numpy as np

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
