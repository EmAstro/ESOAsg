"""
Module to download data from the ESO archive.

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

import os
import sys
import urllib
import numpy as np

from pyvo import dal
from astropy import coordinates
from astropy import units
from astropy import table

from ESOAsg import msgs
from ESOAsg import default
from ESOAsg import ancillary

from IPython import embed

def download(dp_id):
    """
    Parameters
    ----------

    Returns
    -------
    """
    for file_name in dp_id:
        print(file_name)
        # Given a dp_id of a public file, the link to download it is constructed as follows:
        download_url = "http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{}&eso_download=file".format(str(file_name.decode("utf-8")))
        #Files are downloaded in a per-position directory structure.
        # All files matching a given position are store under the directory whose name is the composition of the coordinates (underscore separated).
        # retrieving:
        # last_filepath = file_path
        urllib.request.urlretrieve(download_url, filename=str(file_name.decode("utf-8"))+'.fits')

    msgs.info('Downloading data')
    
    # res['symlink']=''

def query_from_radec(position,
                     instrument=None,
                     maxrec=default.get_value('maxrec')):
    """
    Parameters
    ----------
    position : astropy SkyCoord object
        Coordinates of the sky you want to query in the format
        of an astropy SkyCoord object. Note that at the moment
        it works for one target at the time. For further detail
        see here:
        https://docs.astropy.org/en/stable/coordinates/

    instrument : str
        Name of the instrument to query. Default is None, i.e.
        it will collect all the data

    Returns
    -------
    result_from_query

    """

    # Define TAP SERVICE
    tapobs = dal.tap.TAPService(default.get_value('eso_tap_obs'))
    msgs.info('Querying the ESO TAP service at:')
    msgs.info('{}'.format(str(default.get_value('eso_tap_obs'))))

    # ToDo EMA
    # This should be more flexible and take lists/arrays as input
    if len(position.ra.degree)>1:
        msgs.warning('The position should refer to a single pointing')
        msgs.warning('only the first location is taken into account')
        msgs.working('multi pointing will be available in a future release')
    # Converting from array to number
    RA, Dec = np.float32(position.ra.degree[0]), np.float32(position.dec.degree[0])

    # Define query
    query = """SELECT target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max,
            em_min, dataproduct_type, instrument_name, abmaglim, proposal_id 
            FROM ivoa.ObsCore WHERE obs_release_date < getdate() AND CONTAINS(POINT('',{},{}),
            s_region)=1""".format(RA, Dec)
    msgs.info('The query is:')
    msgs.info('{}'.format(str(query)))

    # Obtaining query results
    result_from_query = tapobs.search(query=query, maxrec=maxrec)

    msgs.info('A total of {} entries has been retrieved'.format(len(result_from_query)))

    return result_from_query
