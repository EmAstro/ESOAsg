from ESOAsg import msgs
import pkg_resources
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from IPython import embed


__all__ = ['Filter']


def _read_filters_file() :
    r"""Read the notes of the table containing the information on the filters currently stored in the filters module
    and store them in an astropy table.

    Returns:
        `astropy.Table`_: Table containing the information on the filters
    """
    all_filters_file = pkg_resources.resource_filename('ESOAsg', 'filters/data/all_filters.txt')
    all_filters_table = ascii.read(all_filters_file, format='commented_header', delimiter=r'\s')
    return all_filters_table


def _read_filter_transmission(file_name, n_filters):
    filter_transmission_file = pkg_resources.resource_filename('ESOAsg', 'filters/data/{}'.format(file_name))
    filter_transmission_table = ascii.read(filter_transmission_file, delimiter=r'\s')
    wavelength = filter_transmission_table['col1'].data
    if n_filters == 1:
        transmission = filter_transmission_table['col2'].data
    elif n_filters == 2:
        transmission = [filter_transmission_table['col2'].data, filter_transmission_table['col3'].data]
    else:
        msgs.error('Number of transmission curves per filter name not implemented')
    return wavelength, transmission


def _read_filters(instrument, mode, filter_name):
    """

    Args:
        instrument:
        mode:
        filter_name:

    Returns:

    """
    filter_table = _read_filters_file()
    mask_instrument = (filter_table['instrument_name'] == instrument)
    mask_mode = (filter_table['instrument_mode'] == mode)
    mask_filter_name = (filter_table['filter_name'] == filter_name)
    mask_combined = mask_instrument & mask_mode & mask_filter_name
    selected_filter_table = filter_table[mask_combined]
    n_filters = selected_filter_table['n_filters'].data[0]
    extracted_file_name = selected_filter_table['file_name'].data[0]
    wavelength, transmission = _read_filter_transmission(extracted_file_name, n_filters)
    return wavelength, transmission


def _read_transmission_curve(instrument, mode, filter_name):
    r"""Given instrument, mode, and filter_name it returns, the transmission curve of the given filter

    Args:
        instrument:
        mode:
        filter_name:

    Returns:

    """
    if instrument is None:
        msgs.warning('Instrument name not defined')
        return None, None
    if mode is None:
        msgs.warning('Instrument mode not defined')
        return None, None
    if filter_name is None:
        msgs.warning('Filter not defined')
        return None, None
    return _read_filters(instrument, mode, filter_name)


def _get_franction_of_transmission(wavelength, transmission, fraction):
    r"""get wavelenght min max given a transmission curve

    Args:
        wavelength:
        transmission:
        fraction:

    Returns:

    """
    peak = np.nanmax(transmission)
    filter_fraction_of_max = (transmission > fraction * peak)
    wavelength_min = np.nanmin(wavelength[filter_fraction_of_max])
    wavelength_max = np.nanmax(wavelength[filter_fraction_of_max])
    return wavelength_min, wavelength_max


class Filter:
    r"""Class that stores the information on a filter
    
    Attributes:
        instrument (str):
        mode (str):
        filter_name (str):
    """

    def __init__(self, instrument=None, mode=None, filter_name=None):
        self.instrument = instrument
        self.mode = mode
        self.filter_name = filter_name

    def get_transmission(self):
        r"""Given a filter it returns the transmission curve

        """
        return _read_transmission_curve(self.instrument, self.mode, self.filter_name)

    def plot(self, show_min_max=True):
        r"""Quick check of the transmission curve of a spectra

        Returns:

        """
        wavelength, transmission = self.get_transmission()
        plt.figure()
        # ToDo
        # this needs to be made more general
        if isinstance(transmission, list):
            plt.plot(wavelength, transmission[0])
            plt.plot(wavelength, transmission[1])
        else:
            plt.plot(wavelength, transmission)
        if show_min_max:
            wavelength_min, wavelength_max = self.get_wavelength_min_max()
            if isinstance(wavelength_min, list):
                for i in np.arange(0,2):
                    plt.axvline(x=wavelength_min[i])
                    plt.axvline(x=wavelength_max[i])
            else:
                plt.axvline(x=wavelength_min)
                plt.axvline(x=wavelength_max)
        return
    
    def get_wavelength_min_max(self, fraction=0.1):
        r"""Return wavelength min max give a filter
        """
        wavelength, transmission = self.get_transmission()
        if isinstance(transmission, list):
            wavelength_min0, wavelength_max0 = _get_franction_of_transmission(wavelength, transmission[0],
                                                                              fraction)
            wavelength_min1, wavelength_max1 = _get_franction_of_transmission(wavelength, transmission[1],
                                                                              fraction)
            return [wavelength_min0, wavelength_min1], [wavelength_max0, wavelength_max1]
        else:
            wavelength_min, wavelength_max = _get_franction_of_transmission(wavelength, transmission,
                                                                            fraction)
            return wavelength_min, wavelength_max