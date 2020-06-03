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

from pyvo import dal
from astropy.coordinates import ICRS
import requests
import webbrowser


from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.ancillary import checks


def _define_tap_service(verbose=False):
    """Load tap service from defaults

    Returns:
        tapcat (`TAP service`)
    """
    if verbose:
        msgs.info('Querying the ESO TAP service at:')
        msgs.info('{}'.format(str(default.get_value('eso_tap_cat'))))
    tapcat = dal.tap.TAPService(default.get_value('eso_tap_cat'))
    return tapcat


def _run_query(query, verbose=False, maxrec=default.get_value('maxrec')):
    tapcat = _define_tap_service()
    if verbose:
        msgs.info('The query is:')
        msgs.info('{}'.format(str(query)))
    # Obtaining query results
    result_from_query = tapcat.search(query=query, maxrec=maxrec)
    return result_from_query


def all_catalogues(verbose=False):
    r"""Load a list with all ESO catalogues

    For further information check `https://www.eso.org/qi/`

    Returns:
        all_catalogues (`list`):
            List of the names of all catalogues currently present at ESO.
    """
    query = """SELECT
                    table_name
               FROM 
                    TAP_SCHEMA.tables 
               WHERE 
                    schema_name='safcat'
             """
    # Obtaining query results
    result_from_query = _run_query(query, verbose=verbose)
    all_catalogues = []
    for catalogue_name in result_from_query['table_name']:
        all_catalogues.append(checks.from_bytes_to_string(catalogue_name))
    all_catalogues.sort()
    return all_catalogues


def query_catalogue(catalogue_name, maxrec=default.get_value('maxrec')):
    r"""Query the ESO tap_cat service (link defined in `ESOAsg\default.txt`) for a specific catalogue.
    
    Args:
        catalogue_name (`str`):
            Catalogue to be queried. To check the full list of catalogues run `all_catalogues()`

    Returns:
        catalogue
    """

    # Check for presence of `catalogue_name` on the ESO archive
    eso_catalogues = all_catalogues(verbose=False)
    if catalogue_name not in eso_catalogues:
        msgs.error('Catalogue: {} not recognized'.format(catalogue_name))

    # query
    query = """SELECT 
                    * 
               FROM 
                    {}
            """.format(catalogue_name)

    # Obtaining query results
    result_from_query = _run_query(query, maxrec=maxrec, verbose=True).to_table()
    return result_from_query
