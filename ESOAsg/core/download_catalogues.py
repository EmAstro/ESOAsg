"""
Module to download data from the ESO catalogues.

------------------- ESO DATA ACCESS POLICY ---------------------
                                                                
The downloaded data are subject to the ESO Data Access Policy   
available at:                                                   
http://archive.eso.org/cms/eso-data-access-policy.html          
                                                                
In particular, you are requested to acknowledge the usage of    
the ESO archive and of the ESO data; please refer to the        
Acknowledgement policies section.                               
                                                                
If you plan to redistribute the downloaded data, please refer   
to the "Requirements for third parties distributing ESO data"   
section.                                                        
                                                                
----------------------------------------------------------------
"""

# import os
# import sys
import urllib
import numpy as np
import os

from astropy.table import MaskedColumn
from pyvo import dal
from astropy.coordinates import ICRS
import requests
import webbrowser


from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.ancillary import checks


def _define_tap_service(verbose=False):
    r"""Load tap service from defaults

    The TAP service for the catalogue is defined in `ESOAsg\default.txt` as `eso_tap_cat`

    Args:
        verbose (`bool`):
            if set to `True` additional info will be displayed

    Returns:
        tapcat (`pyvo.dal.tap.TAPService`)
            TAP service that will be used for the queries
    """
    if verbose:
        msgs.info('Querying the ESO TAP service at:')
        msgs.info('{}'.format(str(default.get_value('eso_tap_cat'))))
    tapcat = dal.tap.TAPService(default.get_value('eso_tap_cat'))
    return tapcat


def _run_query(query, verbose=False, maxrec=default.get_value('maxrec')):
    r"""Run tap query and return result as a table

    Args:
        query (`str`):
            Query to be run
        verbose (`bool`):
            if set to `True` additional info will be displayed
        maxrec (`int`):
            Define the maximum number of entries that a single query can return

    Returns:
        result_from_query (`astropy.Table`):
            Result from the query to the TAP service
    """
    # Load tap service
    tapcat = _define_tap_service(verbose=False)
    if verbose:
        msgs.info('The query is:')
        msgs.info('{}'.format(str(query)))
    # Obtaining query results and convert it to an astropy table
    result_from_query = tapcat.search(query=query, maxrec=maxrec).to_table()
    return result_from_query


def _is_table_at_eso(table_name):
    r"""Check if a given table is present at ESO

    Args:
        table_name (`str`):
            Table to be tested.
    Returns:
        is_at_eso (`bool`):
            `True` if the table is present in tapcat. `False` and warning raised otherwise.
    """
    is_at_eso = True
    # Check for presence of `table_name` on the ESO archive
    eso_catalogues = all_catalogues(verbose=False)['table_name'].data.data.tolist()
    if table_name not in eso_catalogues:
        msgs.warning('Catalogue: {} not recognized. Possible values are:\n{}'.format(table_name, eso_catalogues))
        is_at_eso = False
    return is_at_eso


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


def all_catalogues(verbose=False):
    r"""Load a list with all ESO catalogues

    For further information check `https://www.eso.org/qi/`

    Returns:
        all_catalogues_table (`astropy.table`):
            `astropy.table` containing: `collection`, `table_name`, `title`, `number_rows`, `version`, `acknowledgment`
            of all catalogues currently present at ESO. In addition the column `last_version` is added. This is an
            attempt to remove obsolete catalogues based on the version number and the title of the catalogue.
    """
    query = """SELECT
                    collection, table_name, title, number_rows, version, acknowledgment
               FROM 
                    TAP_SCHEMA.tables 
               WHERE 
                    schema_name='safcat'
             """
    # Obtaining query results
    all_catalogues_table = _run_query(query, verbose=verbose)
    # Sorting
    all_catalogues_table.sort(['collection', 'table_name', 'version'])
    # Checking for obsolete
    unique_titles = np.unique(all_catalogues_table['title'].data).tolist()
    last_version = np.zeros_like(all_catalogues_table['version'].data, dtype=bool)
    for unique_title in unique_titles:
        most_recent_version = np.nanmax(all_catalogues_table['version'].data[(all_catalogues_table[
                                                                                  'title'].data == unique_title)])
        last_version[(all_catalogues_table['title'].data == unique_title) &
                     (all_catalogues_table['version'].data == most_recent_version)] = True
    all_catalogues_table.add_column(MaskedColumn(data=last_version, name='last_version', dtype=bool,
                                                 description='True if this is the latest version of the catalog'))

    all_catalogues_table['table_name'].data.data[:] = checks.from_bytes_to_string(all_catalogues_table[
                                                                                      'table_name'].data.data)
    return all_catalogues_table


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
    query = """SELECT
                    column_name, datatype, description, table_name, unit 
               FROM
                    columns
               WHERE
                    table_name='{}'
            """.format(table_name)
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
    query = """SELECT 
                    {} 
               FROM 
                    {}
            """.format(good_columns_string, table_name)

    # Obtaining query results
    result_from_query = _run_query(query, maxrec=maxrec, verbose=True)
    return result_from_query
