""" Module that performs some lists massaging to work well with arg.parser
"""


import numpy as np
from astropy.table import Column, MaskedColumn
import astropy.units as u

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
            input list of strings that will be checked (usually is coming from `parse_arguments()` in a macro).
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
                list_of_int = [list_of_int[0]] * length
            else:
                msgs.error('List length: {} not matching {}'.format(len(list_of_int), length))
    return list_of_int


def from_element_to_list(element, element_type=str):
    r"""Given an element it returns a list containing the element

    It also checks all the elements in the list have the same type defined by `element_type`

    Args:
        element (any): element that will be put in the list
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
    elif isinstance(element, np.ndarray):
        element_list = element.tolist()
        for element_in_list in element_list:
            assert isinstance(element_in_list, element_type), r'{} must be a {}'.format(element_in_list, element_type)
        return element_list
    elif isinstance(element, MaskedColumn):
        element_list = element.data.data.tolist()
        for element_in_list in element_list:
            assert isinstance(element_in_list, element_type), r'{} must be a {}'.format(element_in_list, element_type)
        return element_list
    elif isinstance(element, element_type):
        return [element]
    else:
        msgs.error('Not valid type for: {}'.format(element))
    return


def from_element_to_list_of_quantities(element, unit=None):
    r"""Convert an input into a list of `astropy.units.Quantity`

    Args:
        element (int, float, np.ndarray, Quantity object, list): element that will be put in the list
        unit (UnitBase instance): An object that represents the unit to be associated with the input value

    Returns:
        list: list of quantities in the format `element`*`unit`

    """
    assert isinstance(unit, u.UnitBase), r'{} not a valid astropy units'.format(unit)
    if isinstance(element, int):
        return [float(element)*unit]
    elif isinstance(element, float):
        return [element*unit]
    elif isinstance(element, np.ndarray) and not isinstance(element, u.Quantity):
        element_list_clean = []
        element_list = np.copy(element).tolist()
        for element_in_list in element_list:
            element_list_clean.append(element_in_list*unit)
        return element_list_clean
    elif isinstance(element, u.Quantity):
        element_list_clean = []
        element_converted = np.copy(element).to(unit)
        for element_in_list in np.nditer(element_converted):
            print(element_in_list)
            element_list_clean.append(element_in_list*unit)
        return element_list_clean
    elif isinstance(element, list):
        element_list_clean = []
        for element_in_list in element:
            if isinstance(element_in_list, u.Quantity):
                element_list_clean.append(element_in_list.to(unit))
            else:
                element_list_clean.append(element_in_list*unit)
        return element_list_clean
    else:
        msgs.error('The input cannot be converted into a list of quantities')
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
        input_in_bytes (any): input in bytes

    Returns:
        any: output converted into a string

    """
    # condition for a single entry as bytes
    if isinstance(input_in_bytes, bytes):
        output_as_str = str(input_in_bytes.decode("utf-8"))
    # condition for a np.array containing bytes
    elif isinstance(input_in_bytes, np.ndarray) and len(input_in_bytes.shape) == 1:
        output_as_str = np.empty_like(input_in_bytes)
        for idx in np.arange(0, np.size(output_as_str)):
            if isinstance(input_in_bytes[idx], bytes):
                output_as_str[idx] = str(input_in_bytes[idx].decode("utf-8"))
            else:
                output_as_str[idx] = input_in_bytes[idx]
    # condition for a list containing bytes
    elif isinstance(input_in_bytes, list):
        output_as_str = []
        for element_in_bytes in input_in_bytes:
            if isinstance(element_in_bytes, bytes):
                output_as_str.append(str(element_in_bytes.decode("utf-8")))
            else:
                output_as_str.append(element_in_bytes)
    # return copy of the entry
    else:
        output_as_str = input_in_bytes.copy()
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
