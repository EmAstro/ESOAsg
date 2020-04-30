"""
Module to perform some quick astronomical calculations
"""

import numpy as np
import requests
from astropy import units as u
from astropy import constants

from astroquery.mast import Mast
from astroquery.mast import Observations

from ESOAsg import msgs
from ESOAsg.ancillary import checks


def download_kepler_data(target_name):
    r"""Run a query to MAST to obtain the kepler data of an object. Date are saved in the directories:
    `./mastDownload/Kepler/<target_name>_*/*.fits`

    Args:
        target_name (`str`):
            Name of the target: e.g., 'kplr011446443'

    Returns:
        The code returns `True` if data are retrieved and `False` otherwise.

    """
    keplerObs = Observations.query_criteria(target_name=target_name, obs_collection='Kepler')
    if len(keplerObs) == 0:
        msgs.warning('No Kepler data for target {}'.format(target_name))
        return False
    for idx in np.arange(0, len(keplerObs)):
        keplerProds = Observations.get_product_list(keplerObs[idx])
        keplerFitsProds = Observations.filter_products(keplerProds, extension='fits', mrp_only=False)
        Observations.download_products(keplerFitsProds, mrp_only=False, cache=False)
    return True


def download_gw_bayestar(superevent_name):
    r"""Download the bayestar.fits.gz of a GW superevent

    This is simply grabbing the file from:
    `https://gracedb.ligo.org/superevents/<superevent_name>/files/`

    Args:
        superevent_name (`str`):
            Name of the superevent: e.g., 'S191205ah'

    Returns:
        The code returns `True` if data are retrieved and `False` otherwise.

    """
    # https://gracedb.ligo.org/superevents/S191205ah/files/
    gw_url = 'https://gracedb.ligo.org/superevents/S191205ah/files/'+superevent_name+'/files/bayestar.fits.gz'
    checks.connection_to_website(gw_url)
    r = requests.get(gw_url, allow_redirects=True)
    open(superevent_name+'_bayestar.fits.gz', 'wb').write(r.content)
    return True

# ToDo
def abmaglim(rms, seeing_fwhm, exptime=None, zero_point=None, sigma=5.):
    r"""Calculate the N-sigma magnituide limit.

    The code behaves in different ways depending if zero_point is set to None or not.
    If zero_points is None:
        rms needs to have units and exptime is not considered.
        abmaglim = -2.5 * Log10( sigma * rms * PI*(.5*seeing_fwhm)**2. / (3631.*u.jansky) )
    else:
        abmaglim is calculated as:
        abamglim = -2.5 * Log10( sigma * rms * PI*(.5*seeing_fwhm)**2. / exptime ) + zero_point

    Args:
        rms (`float`):
           Calculated rms on the image
        exptime (`float`):
            Exposure time in seconds
        zero_point (`float`):
            AB zero point of the image. If not set, the code will take the RMS with the associated
            astropy.units.
        seeing_fwhm (`float`):
            FWHM of the seeing in pixels
        sigma (`float`):
           Number of sigma to consider a detection significant

    Returns:
        abmaglim (`float`):
            AB magnitude limit of for an image
    """
    msgs.work('Calculating {}-sigma AB mag limit'.format(sigma))

    if zero_point is None:
        if exptime is not None:
            msgs.error('`exptime` needs to be `None` if `zero_point` is `None`')
        rms_in_jansky = rms.to(u.jansky, equivalencies=u.spectral())
        abmaglim = -2.5*np.log10(sigma*rms_in_jansky*np.pi*np.power(seeing_fwhm/2.,2.)/(3631.*u.jansky))
    else:
        abmaglim = -2.5 * np.log10(sigma*rms*np.pi*np.power(seeing_fwhm/2.,2.)/exptime) + zero_point
    return abmaglim

# ToDo
def fnu2flambda(flambda, wavelength):
    return np.power(wavelength, 2.)*flambda/constants.c

# ToDo
def flambda2fnu(fnu, wavelength):
    return constants.c*fnu/np.power(wavelength, 2.)
