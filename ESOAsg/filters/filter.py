from ESOAsg import msgs
import pkg_resources
from astropy.io import ascii

__all__ = ['Filter']


def _read_filters_file() -> dict:
    r"""Read the notes of the table containing the information on the filters currently stored in the filters module
    and store them in an astropy table.

    Returns:
        `astropy.Table`_: Table containing the information on the filters
    """
    all_filters_file = pkg_resources.resource_filename('ESOAsg', 'filters/data/all_filters.txt')
    all_filters_table = ascii.read(all_filters_file, format='commented_header', delimiter=r'\s')
    return all_filters_table


class Filter:
    r"""Class that stores the information on a filter
    
    Attributes:
        instrument (str):
        mode (str):
        filter_name (str):
        wavelmin (float):
        wavelmax (float):
    """

    def __init__(self, instrument=None, mode=None, filter_name=None):
        self.instrument = instrument
        self.mode = mode
        self.filter_name = filter_name
