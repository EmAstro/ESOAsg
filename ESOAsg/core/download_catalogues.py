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
        query (`str, np.str`):
            Query to be run
        verbose (`bool`):
            if set to `True` additional info will be displayed
        maxrec (`int, numpy.int`):
            Define the maximum number of entries that a single query can return

    Returns:
        result_from_query (`astropy.table`):
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
                                                                                  'title'].data==unique_title)])
        last_version[(all_catalogues_table['title'].data==unique_title) &
                     (all_catalogues_table['version'].data==most_recent_version)] = True
    all_catalogues_table.add_column(MaskedColumn(data=last_version, name='last_version', dtype=bool,
                                                 description='True if this is the latest version of the catalog'))
    return all_catalogues_table


def rows_in_catalogue(catalogue_name):
    r"""

    Args:
        catalogue_name (`str`):
            Catalogue to be queried. To check the full list of catalogues run `all_catalogues()`

    Returns:

    """


def query_catalogue(table_name, maxrec=default.get_value('maxrec')):
    r"""Query the ESO tap_cat service (link defined in `ESOAsg\default.txt`) for a specific catalogue.
    
    Args:
        table_name (`str`):
            Table to be queried. To check the full list of catalogues run `all_catalogues()`
        maxrec (`int, numpy.int`, `None`):
            Define the maximum number of entries that a single query can return. If set to `None` the
            entire catalogues is queried (this may cause size problems). The default values is set in
            `ESOAsg\default.txt` as `max_rec`

    Returns:
        catalogue
    """

    # Check for presence of `table_name` on the ESO archive
    eso_catalogues = all_catalogues(verbose=False)
    test_table_names = [checks.from_bytes_to_string(test_table_name) for test_table_name in eso_catalogues[
        'table_name'].data.tolist()]
    if table_name not in test_table_names:
        msgs.error('Catalogue: {} not recognized. Possible values are:\n{}'.format(table_name, test_table_names))

    # query
    query = """SELECT 
                    * 
               FROM 
                    {}
            """.format(table_name)

    # Obtaining query results
    result_from_query = _run_query(query, maxrec=maxrec, verbose=True)
    return result_from_query
