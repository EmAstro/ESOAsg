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
        print(str(file_name.decode("utf-8"))
        # Given a dp_id of a public file, the link to download it is constructed as follows:
        download_url = "http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{}&eso_download=file".format(str(file_name.decode("utf-8")))
        print(download_url)
        #Files are downloaded in a per-position directory structure.
        # All files matching a given position are store under the directory whose name is the composition of the coordinates (underscore separated).
        # retrieving:
        # last_filepath = file_path
        urllib.request.urlretrieve(download_url)

    msgs.info('Downloading data')
    
    # res['symlink']=''

def query_from_radec(position,
                     maxrec=default.get_value('maxrec')):
    """
    Parameters
    ----------
    position : astropy SkyCoord object
        Coordinates of the sky you wnat to query in the format
        of an astropy SkyCoord object. For further detail see
        here: https://docs.astropy.org/en/stable/coordinates/

    Returns
    -------
    result_from_query

    """
    # Define TAP SERVICE
    tapobs = dal.tap.TAPService(default.get_value('eso_tap_obs'))
    msgs.info('Querying the ESO TAP service at:')
    msgs.info('{}'.format(str(default.get_value('eso_tap_obs'))))

    # converting from array to number
    RA, Dec = np.float32(position.ra.degree[0]), np.float32(position.dec.degree[0])

    # Define query
    query = """SELECT dp_id, s_ra, s_dec, em_min*1E9 min_wave_nm, em_max*1E9 max_wave_nm, t_min, abmaglim, proposal_id, access_estsize FROM ivoa.ObsCore WHERE obs_release_date < getdate() AND obs_collection='MUSE' AND CONTAINS(POINT('',{},{}), s_region)=1""".format(RA, Dec)

    print(query)

    result_from_query = tapobs.search(query=query, maxrec=maxrec)

    return result_from_query
