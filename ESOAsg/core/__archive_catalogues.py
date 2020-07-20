
# import os
# import sys
# import urllib
import numpy as np
# import os

from astropy.table import MaskedColumn
from pyvo import dal
# from astropy.coordinates import ICRS
# import requests
# import webbrowser


from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.ancillary import checks





def _are_columns_in_table(column_list, table_name):
    r"""Check if a given table is present at ESO

    Args:
        column_list (`list`):
            list of the columns that will be checked
        table_name (`str`):
            Table to be tested.
    Returns:
        are_in_table (`list`):
            `True` if the column is present in the table. `False` and warning raised otherwise.
    """
    are_in_table = []
    first_error = True
    eso_column_list = columns_in_catalogue(table_name, verbose=False)['column_name'].data.data.tolist()
    for column_test in column_list:
        if column_test in eso_column_list:
            are_in_table.append(True)
        else:
            are_in_table.append(False)
            if first_error:
                msgs.warning('Column: {} not recognized. Possible values are:\n{}'.format(column_test, eso_column_list))
                first_error = False
            else:
                msgs.warning('Column: {} not recognized'.format(column_test))
    return are_in_table

def columns_in_catalogue(table_name, verbose=False):
    r"""Return the columns present in a table

    Args:
        table_name (`str`):
            Table to be queried. To check the full list of catalogues run `all_catalogues()`
        verbose (`bool`):
            if set to `True` additional info will be displayed

    Returns:
        columns_table (`astropy.table`):
            List of all columns present in a table. Information are stored in `column_name`, `datatype`, `description`,
            and `unit`
    """
    # Check for the table
    if not _is_table_at_eso(table_name):
        msgs.error('{} is not a valid table'.format(table_name))
    query = '''
            SELECT
                column_name, datatype, description, table_name, unit 
            FROM
                columns
            WHERE
                table_name='{}'
            '''.format(table_name)
    # Obtaining query results
    columns_table = _run_query(query, maxrec=None, verbose=verbose)
    columns_table['column_name'].data.data[:] = checks.from_bytes_to_string(columns_table['column_name'].data.data)
    columns_table.remove_column('table_name')
    return columns_table


def query_catalogue(table_name, which_columns=None, maxrec=default.get_value('maxrec')):
    r"""Query the ESO tap_cat service (link defined in `ESOAsg\default.txt`) for a specific catalogue.
    
    Args:
        table_name (`str`):
            Table to be queried. To check the full list of catalogues run `all_catalogues()`
        which_columns (`list`):
            List of the columns that you want to download. The full list of the columns in a table can be found
            running `columns_in_catalogue(table_name)`
        maxrec (`int`):
            Define the maximum number of entries that a single query can return. If set to `None` the
            value is set by the limit of the service. The default values is set in`ESOAsg\default.txt` as `max_rec`

    Returns:
        catalogue
    """
    # Check for the table
    if not _is_table_at_eso(table_name):
        msgs.error('{} is not a valid table'.format(table_name))

    # Check for columns and select the one present in the table
    if which_columns is not None:
        test_columns = _are_columns_in_table(which_columns, table_name)
        good_columns = [which_column for which_column, test_column in zip(which_columns, test_columns) if test_column]
        bad_columns = [which_column for which_column, test_column in zip(which_columns, test_columns)
                       if not test_column]
        if len(bad_columns) > 0:
            msgs.warning('The following columns will be excluded from the query: {}'.format(bad_columns))
        if len(good_columns) == 0:
            good_columns = ['*']
            msgs.warning('No valid column selected. All the columns will be queried')
    else:
        good_columns = ['*']

    # Create string of all columns for the query
    good_columns_string = str(' ')
    for good_column in good_columns:
        good_columns_string = good_columns_string + str(good_column) + ', '
    good_columns_string = good_columns_string.strip()[:-1]
    # query
    query = '''
            SELECT 
                {} 
            FROM 
                {}
            '''.format(good_columns_string, table_name)

    # Obtaining query results
    result_from_query = _run_query(query, maxrec=maxrec, verbose=True)
    return result_from_query


def query_from_radec(table_names, position, radius=None, maxrec=default.get_value('maxrec')):
    r"""Query the ESO tap_cat service (link defined in `ESOAsg\default.txt`) for a specific a specific location.

    The `position` needs to be given as an `astropy.coordinates.SkyCoord` object.

    Args:
        table_names (`list`):
            list of the tables to be queried. To check the full list of catalogues run `all_catalogues()`
        position (`astropy.coordinates.SkyCoord`):
            Coordinates of the sky you want to query in the format of an `astropy.coordinates.SkyCoord` object. Note
            that at the moment it works for one target at the time. For further detail see here:
            `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_
        radius (`float`):
            Search radius you want to query in arcseconds. Note that in case `None` is given, the query will be
            performed with the `CONTAINS(POINT('',RA,Dec), s_region)` clause instead of the
            `CONTAINS(s_region,CIRCLE('',RA,Dec,radius/3600.))` one. See here for further examples:
            `tap obs examples <http://archive.eso.org/tap_obs/examples>`_
        maxrec (`numpy.int`):
            Define the maximum number of file that a single query can return from the ESO archive. You probably never
            need this. By default is set by the `default.txt` file.

    Returns:
        result_from_query (`astropy.Table`):
            Result from the query to the TAP service
    """