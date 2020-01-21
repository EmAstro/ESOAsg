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


def download():
    """
    Parameters
    ----------

    Returns
    -------
    """
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
    """
    # Define TAP SERVICE
    tapobs = dal.tap.TAPService(default.get_value('eso_tap_obs'))
    msgs.info('Querying the ESO TAP service at:')
    msgs.info('{}'.format(str(default.get_value('eso_tap_obs'))))

    RA, Dec = np.float32(position.ra.degree), np.float32(position.dec.degree)

    # Define query
    query = """SELECT J0100 as target, dp_id, s_ra, s_dec, em_min*1E9 min_wave_nm, em_max*1E9 max_wave_nm, t_min, abmaglim, proposal_id, access_estsize FROM ivoa.ObsCore WHERE obs_release_date < getdate() AND obs_collection='MUSE' AND CONTAINS(POINT('',{},{}), s_region)=1""".format(RA, Dec)

    print(query)

    res = tapobs.search(query=query, maxrec=maxrec)
    tab = res.to_table()
    print(tab)
