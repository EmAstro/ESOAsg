"""
Module to perform some quick astronomical calculations
"""

import numpy as np
from astropy import units as u
from astropy import constants
from ESOAsg import msgs


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

def fnu2flambda(flambda, wavelength):
    return np.power(wavelength, 2.)*flambda/constants.c

def flambda2fnu(fnu, wavelength):
    return constants.c*fnu/np.power(wavelength, 2.)
