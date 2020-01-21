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


from pyvo import dal
from astropy import coordinates
from astropy import units


from ESOAsg import msgs
from ESOAsg import default
from ESOAsg import ancillary


def download():
    # res['symlink']=''
    print(default.get_value('eso_tap_obs'))

def query():
    """
    """
    # Define TAP SERVICE
    tapobs = dal.tap.TAPService(default.get_value('eso_tap_obs'))
    msgs.info('Querying the ESO TAP service at {:.10}'.format(str(default.get_value('eso_tap_obs'))))
