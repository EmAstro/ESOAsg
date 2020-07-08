""" Module that performs some lists massaging to work well with arg.parser
"""

from ESOAsg import msgs
from ESOAsg.ancillary import checks


def make_list_of_fits_files(args_input, length=None):
    r"""Cleaning an input list of fits files

    Args:
        args_input ('list'):
            input list of fits files that will be checked (usually is coming from
            `parse_arguments()` in a macro).
    Returns:
        list_of_fits_files ('list'):
            list containing all the valid fits files given in input
    """
    list_of_fits_files = []
    if not isinstance(args_input, list):
        args_input_files = [args_input]
    else:
        args_input_files = args_input
    for args_input_file in args_input_files:
        if checks.fits_file_is_valid(args_input_file, overwrite=False, verify_fits=False):
            list_of_fits_files.append(args_input_file)
        else:
            msgs.warning('{} excluded because not a valid fits file'.format(args_input_file))
    if len(list_of_fits_files) == 0:
        msgs.error('No valid fits files present')
    return list_of_fits_files


def make_list_of_strings(args_input, length=None):
    r"""Cleaning an input list of strings

    Args:
        args_input (`list`):
            input list of strings that will be checked (usually is coming from
            `parse_arguments()` in a macro).
        length (`int`):
            If set to a value, the code will check that the length of `list_of_string` will match `length`. If not
            the code will further check if `args_input` contains only one element. In this case the output list will
            contain this values `length` times. If this situation does not happen, an error is raised.
    Returns:
        list_of_string (`list`):
            list containing all the valid strings given in input
    """
    if length is not None:
        assert(isinstance(length, int)), '`length` must be an integer'
    list_of_string = []
    if not isinstance(args_input, list):
        args_input_strings = [args_input]
    else:
        args_input_strings = args_input
    for args_input_string in args_input_strings:
        if isinstance(args_input_string, str):
            list_of_string.append(args_input_string)
        else:
            msgs.warning('{} excluded because not a valid string'.format(args_input_string))
    if len(list_of_string) == 0:
        msgs.error('No valid element present in the list')
    if length is not None:
        if len(list_of_string) != length:
            if len(list_of_string) == 1:
                list_of_string = [list_of_string[0]] * length
            else:
                msgs.error('List length: {} not matching {}'.format(len(list_of_string), length))
    return list_of_string


def make_list_of_int(args_input, length=None):
    r"""Cleaning an input list of int

    Args:
        args_input (`list`):
            input list of strings that will be checked (usually is coming from
            `parse_arguments()` in a macro).
        length (`int`):
            If set to a value, the code will check that the length of `list_of_string` will match `length`. If not
            the code will further check if `args_input` contains only one element. In this case the output list will
            contain this values `length` times. If this situation does not happen, an error is raised.
    Returns:
        list_of_int (`list`):
            list containing all the valid int given in input
    """
    if length is not None:
        assert(isinstance(length, int)), '`length` must be an integer'
    list_of_int = []
    if not isinstance(args_input, list):
        args_input_strings = [args_input]
    else:
        args_input_strings = args_input
    for args_input_string in args_input_strings:
        if isinstance(args_input_string, int):
            list_of_int.append(args_input_string)
        else:
            msgs.warning('{} excluded because not a valid string'.format(args_input_string))
    if len(list_of_int) == 0:
        msgs.error('No valid element present in the list')
    if length is not None:
        if len(list_of_int) != length:
            if len(list_of_int) == 1:
                list_of_string = [list_of_int[0]] * length
            else:
                msgs.error('List length: {} not matching {}'.format(len(list_of_int), length))
    return list_of_int
