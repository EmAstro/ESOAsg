r"""Class to work on light curves
"""

import numpy as np
from astropy.io import fits
from astropy.io.fits.column import NUMPY2FITS
from astropy.table import Column, Table

from ESOAsg import msgs
from ESOAsg.core import fitsfiles
from ESOAsg.ancillary import checks


def _return_data_from_column(table, col_name):
    r"""Check that data are present in the `col_name` column of a table.

    Args:
        table (`fits.BinTableHDU`, `fits.TableHDU`):
            table from which you want to extract the data
        col_name (`str`):
            column name from which you want to extract the data.

    Returns:
        column_data (np.array):
            numpy array containing the data. Set to `None` if no data are present
        column_length (np.int):
            length of column_data. Set to `np.int(0)` if no data are present
    """
    if len(table.data[col_name]) > 0:
        column_data = np.array(table.data[col_name])
        column_length = len(column_data)
    else:
        column_data = None
        column_length = np.int(0)
    return column_data, column_length


def _return_attribute_from_column(table, col_name, attribute):
    r"""Returns the value of an attribute from a column of a table

    Args:
        table (`fits.BinTableHDU`, `fits.TableHDU`):
            table from which you want to extract the information.
        col_name (`str`):
            column name from which you want to extract the information.
        attribute (`str`):
            attribute name

    Returns:
        attribute_value:
            value of the selected attribute. Returns `None` if attribute is not present
    """
    if hasattr(table.columns[col_name], attribute):
        attribute_value = getattr(table.columns[col_name], attribute)
    else:
        attribute_value = None
    return attribute_value


def _check_attribute_is_in_light_curve_columns(attribute): # Written by Ema 16.04.2020
    r"""Check if an attribute of a `LightCurves` object is not `header`, `primary_header`, `other`, or `_datatype`

    Args:
        attribute (`str`):
            attribute of the `LightCurves` object you want to check

    Returns:
        is_in_column (`bool`):
            `True` if the `attribute` is not `header`, `primary_header`, `other`, or `_datatype`

    """
    assert isinstance(attribute, (str, np.str)), 'Not valid string for `attribute`'
    is_in_column = not attribute.endswith('header') and attribute is not 'others' and not attribute.startswith('_')
    return is_in_column


def _return_format_column(column):
    r"""Returns format of a column in the fits format

    see https://docs.astropy.org/en/stable/io/fits/usage/table.html : column creation for further details.

    Args:
        column (`astropy.table.Column`):
            a `Column` object

    Returns:
        column_format (`str`):
            a string containing the size and the dtype of of `column`
    """
    assert isinstance(column, Column), r"Not a valid `column`"
    column_size = str(column.size)
    #
    if column.dtype.str.startswith(('<','>','=')):
        column_dtype = column.dtype.str[1:]
    else:
        column_dtype = column.dtype.str
    column_type = NUMPY2FITS[column_dtype]
    column_format = column_size + column_type
    return column_format

def save_into_fits(fits_file_name, primary_header, light_curve_headers, light_curve_names, light_curves,
                   overwrite=True):
    r"""

    Args:
        light_curves:
        light_curve_names:
        light_curve_headers:
        fits_file_name:
        primary_header:
        lightcurves:
        overwrite (`bool`):

    Returns:

    """
    assert isinstance(fits_file_name, (str, np.str)), r'Not a valid name for the fits file'
    if not fits_file_name.endswith('.fits'):
        fits_files_name = fits_file_name + '.fits'

    if not isinstance(light_curves, list):
        assert isinstance(light_curves, LightCurves), r'Not a LightCurves object'
        light_curves_list = [light_curves]
    else:
        for light_curve in light_curves:
            assert isinstance(light_curve, LightCurves), r'Not a LightCurves object'
        light_curves_list = light_curves

    if light_curve_names is None:
        light_curves_names_list = [None] * len(light_curves_list)
    elif not isinstance(light_curve_names, list):
        assert isinstance(light_curve_names, (str, np.str)), r'Not a valid name'
        light_curves_names_list = [light_curve_names]
    else:
        for light_curve_name in light_curve_names:
            assert isinstance(light_curve, (str, np.str)), r'Not a valid name'
        light_curves_names_list = light_curve_name

    if light_curve_headers is None:
        light_curve_headers_list = [None] * len(light_curves_list)
    elif not isinstance(light_curve_headers, list):
        assert isinstance(light_curve_headers, fits.Header), r'Not a valid header'
        light_curve_headers_list = [light_curve_headers]
    else:
        for light_curve_header in light_curve_headers:
            assert isinstance(light_curve_header, fits.Header), r'Not a valid header'
        light_curve_headers_list = light_curves

    if len(light_curves_names_list) != len(light_curves_list):
        msgs.error(r'`light_curve_names` and `light_curves` must have the same length')

    if len(light_curve_headers_list) != len(light_curves_list):
        msgs.error(r'`light_curve_headers` and `light_curves` must have the same length')

    # Primary HDU
    primary_hdu = fits.PrimaryHDU()
    if primary_header is not None:
        assert isinstance(primary_header, fits.Header), r'Not a valid header'
        primary_hdu.header = primary_header

    hdul = fits.HDUList(primary_hdu)

    # saving the lightcurves:
    for light_curve_name, light_curve, light_curve_header in zip(light_curves_names_list, light_curves_list, light_curve_headers_list):
        msgs.info('Saving {}'.format(light_curve_name))
        list_of_columns = []
        for attribute in light_curve.__dict__.keys():
            if _check_attribute_is_in_light_curve_columns(attribute):
                if getattr(light_curve, attribute) is not None:
                    # Name
                    col_name = attribute.upper()
                    # Format
                    col_format = _return_format_column(getattr(light_curve, attribute))
                    # Units
                    # ToDo Units
                    col_units = getattr(light_curve, attribute).unit
                    # Data
                    col_array = np.array(getattr(light_curve, attribute).data)
                    # ToDo Units
                    col = fits.Column(name=col_name, format=col_format, array=[col_array])
                    list_of_columns.append(col)
        if light_curve.others is not None:
            for other_attribute in light_curve.others.colnames:
                # Name
                col_name = other_attribute.upper()
                # Format
                col_format = _return_format_column(light_curve.others[other_attribute])
                # Units
                # ToDo Units
                col_units = light_curve.others[other_attribute].unit
                # Data
                col_array = np.array(light_curve.others[other_attribute].data)
                # ToDo Units
                col = fits.Column(name=col_name, format=col_format, array=[col_array])
                list_of_columns.append(col)
        if light_curve_header is not None:
            hdu = fits.BinTableHDU.from_columns(list_of_columns, header=light_curve_header, nrows=1)
        else:
            hdu = fits.BinTableHDU.from_columns(list_of_columns,nrows=1)
        if light_curve_name is not None:
            hdu.header['EXTNAME'] = light_curve_name
        hdul.append(hdu)
    hdul.writeto(fits_file_name, overwrite=overwrite)


class LightCurves:
    r"""A class used to define and make simple operations on time series

    This allows to perform some basic tasks on a time series and to save it in a format that is ESO Phase3
    compliant.

    Attributes:

    Methods:

    """

    def __init__(self, primary_header=None, header=None, time=None, flux=None, error=None, time_bin=None,
                 background=None, quality=None, others=None):
        r"""Instantiate the class LightCurves


        Each field of the BINTABLE shall be further described in the extension header. Mandatory fields shall be:
        * time
        * flux
        * error
        in that particular order.
        Additional fields may be added.



        """
        self.primary_header = primary_header
        self.header = header
        self.time = time
        self.time_bin = time_bin
        self.flux = flux
        self.error = error
        self.background = background
        self.quality = quality
        if others is not None:
            self.others = others
        else:
            self.others = Table()
        self._datatype = 'LightCurves'

    def load_from_table(self, table, primary_header=None, copy_header=True, where_time='TIME',
                        where_time_bin='TIME_BIN', where_flux='FLUX',
                        where_error='ERROR', where_background='BACKGROUND', where_quality='QUAL'):
        r"""Given a table put it in a LightCurves object

        Args:
            where_quality:
            where_background:
            where_error:
            where_time_bin:
            where_time:
            copy_header:
            primary_header:
            where_flux :
        """
        if checks.table_is_valid(table):
            msgs.work('Reading input table')

        if primary_header is not None:
            if len(primary_header) > 0:
                self.primary_header = primary_header
            else:
                msgs.warning('Empty `primary_header` provided')

        if copy_header:
            if len(table.header) > 0:
                self.header = table.header
            else:
                msgs.warning('No header found in the table')

        if isinstance(table, fits.BinTableHDU):
            self._load_from_BinTableHDU(table, copy_header=copy_header, where_time=where_time,
                                        where_time_bin=where_time_bin, where_flux=where_flux, where_error=where_error,
                                        where_background=where_background, where_quality=where_quality)
        elif isinstance(table, fits.TableHDU):
            # ToDo implement TableHDU case
            msgs.error('To be implemented')
        else:
            msgs.error('Unknown table type')

    def _load_from_BinTableHDU(self, table, copy_header=True, where_time='TIME', where_time_bin='TIME_BIN',
                               where_flux='FLUX', where_error='ERROR', where_background='BACKGROUND',
                               where_quality='QUAL'):
        r"""Given a BinTableHDU put it in a LightCurves object
        """

        # Going through all the columns names
        column_names = [column.name for column in table.columns]
        column_loaded = [False for column in table.columns]

        #  Loading time, time_bin, flux, error, background, qual (if possible)
        for attribute, value in zip(self.__dict__.keys(), self.__dict__.values()):
            if _check_attribute_is_in_light_curve_columns(attribute):
                where_attribute = vars()['where_' + attribute]

                # Loading attributes from table columns
                if where_attribute in column_names:
                    msgs.info('Loading {} from column {}'.format(attribute, where_attribute))

                    # Check that data are present
                    attribute_data, attribute_length = _return_data_from_column(table, where_attribute)
                    if attribute_length == 0:
                        msgs.warning('No data present. {} set to `None`'.format(attribute))

                    # Checks for units and dtype
                    attribute_unit = _return_attribute_from_column(table, where_attribute, 'unit')
                    attribute_dtype = _return_attribute_from_column(table, where_attribute, 'dtype')

                    attribute_column = Column(name=attribute, data=attribute_data, length=attribute_length,
                                              unit=attribute_unit, dtype=attribute_dtype)
                    setattr(self, attribute, attribute_column)
                    column_loaded[column_names.index(where_attribute)] = True

        # Dealing with others
        for other_attribute in column_names:
            if column_loaded[column_names.index(other_attribute)] is False:
                other_attribute_data, other_attribute_length = _return_data_from_column(table, other_attribute)
                other_attribute_unit = _return_attribute_from_column(table, other_attribute, 'unit')
                other_attribute_dtype = _return_attribute_from_column(table, other_attribute, 'dtype')
                other_attribute_column = Column(name=other_attribute, data=other_attribute_data,
                                                length=other_attribute_length, unit=other_attribute_unit,
                                                dtype=other_attribute_dtype)
                self.others.add_column(other_attribute_column)
                column_loaded[column_names.index(other_attribute)] = True

    def check(self):
        r"""Checks that a LightCurves objects is in a format compatible with the ESO standard
        """
        good_light_curve = True

        # Checks that time, flux, and error are inside the LightCurve object
        if self.time is None:
            msgs.warning('`time` is not defined')
            good_light_curve = False
        if self.flux is None:
            msgs.warning('`flux` is not defined')
            good_light_curve = False
        if self.error is None:
            msgs.warning('`error` is not defined')
            good_light_curve = False

        # Check that all columns have the same length
        test_length = self.time.size
        for attribute in self.__dict__.keys():
            if _check_attribute_is_in_light_curve_columns(attribute):
                if getattr(self, attribute) is not None:
                    if getattr(self, attribute).size != test_length:
                        msgs.warning('Inconsistent length in {}'.format(attribute))
                        good_light_curve = False
        for other_attribute in self.others.colnames:
            if self.others[other_attribute].size != test_length:
                msgs.warning('Inconsistent length in `others[{}]`'.format(attribute))
                good_light_curve = False

        # Check that time is monotonic and always positive
        if np.min(self.time.data) < 0.:
            msgs.warning('`time` must be positive')
            good_light_curve = False
        delta_time = np.diff(self.time.data)
        if np.any(delta_time <= 0):
            msgs.warning('`time` should increase monotonically')
            good_light_curve = False

        # Check that there are no +/-inf
        for attribute in self.__dict__.keys():
            if _check_attribute_is_in_light_curve_columns(attribute):
                if getattr(self, attribute) is not None:
                    if np.any(np.isinf(getattr(self, attribute).data)):
                        msgs.warning('{} contains +/-inf values'.format(attribute))
                        good_light_curve = False
        for other_attribute in self.others.colnames:
            if np.any(np.isinf(self.others[other_attribute].data)):
                msgs.warning('`others[{}]` contains +/-inf values'.format(other_attribute))
                good_light_curve = False

        return good_light_curve

    def to_fits(self, fits_file_name, light_curve_name='LIGHTCURVE', overwrite=True):
        if not self.check():
            msgs.error('the LightCurve object does not respect the requirements from ESO')
        save_into_fits(fits_file_name, self.primary_header, self.header, light_curve_name, self, overwrite=overwrite)
