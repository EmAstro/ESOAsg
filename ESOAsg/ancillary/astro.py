"""
Module to perform some quick astronomical calculations

The GW part is taken from a python script written by Martino Romaniello
http://www.eso.org/~mromanie/martino/Home.html
"""

import numpy as np
import requests
import json

import matplotlib
from matplotlib import pyplot as plt
STARTING_MATPLOTLIB_BACKEND = matplotlib.rcParams["backend"]

from astropy import units as u
from astropy import constants

import importlib

from os import path
from os import remove
from astroquery.mast import Mast
from astroquery.mast import Observations
from ligo.skymap.tool.ligo_skymap_contour import main as ligo_skymap_contour
matplotlib.rcParams["backend"] = STARTING_MATPLOTLIB_BACKEND

import healpy
from astropy import wcs


from ESOAsg import msgs
from ESOAsg.ancillary import checks
from ESOAsg.ancillary import polygons # ToDo to be removed
from ESOAsg.core import fitsfiles


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


def download_gw_bayestar(superevent_name, file_name='bayestar.fits.gz'):
    r"""Download the bayestar.fits.gz of a GW superevent

    This is simply checking the existence of:
    `https://gracedb.ligo.org/superevents/<superevent_name>/files/`
    and downloading the file:
    `https://gracedb.ligo.org/apiweb/superevents/<superevent_name>/iles/bayestar.fits.gz`

    The downloaded file will be renamed as `superevent_name`+`file_name`

    Args:
        superevent_name (`str`):
            Name of the superevent: e.g., 'S191205ah'
        file_name (`str`):
            Name of the file to download. Default is `bayestar.fits.gz`

    Returns:
        The code returns `True` if data are retrieved and `False` otherwise.
    """

    # Some checks
    assert isinstance(superevent_name, (str, np.str)), '{} is not a valid string'.format(superevent_name)
    assert isinstance(file_name, (str, np.str)), '{} is not a valid string'.format(file_name)
    # Checks if superevent exists
    gw_files_url = 'https://gracedb.ligo.org/superevents/'+superevent_name+'/files/'
    checks.connection_to_website(gw_files_url)

    # download the and save the superevent file
    gw_bayestar_url = 'https://gracedb.ligo.org/apiweb/superevents/'+superevent_name+'/files/'+file_name
    r = requests.get(gw_bayestar_url, allow_redirects=True)
    if r.status_code == 404:
        msgs.warning('Failed to access to: {}'.format(gw_bayestar_url))
        return False
    open(superevent_name+'_bayestar.fits.gz', 'wb').write(r.content)

    # check that the file actually arrived on disk
    if path.isfile(superevent_name+'_bayestar.fits.gz'):
        _test_header = fitsfiles.header_from_fits_file(superevent_name + '_' + file_name)
        if len(_test_header) > 0:
            msgs.info('File {}_{} successfully downloaded'.format(superevent_name, file_name))
        else:
            msgs.warning('Not a valid fits file in: {}'.format(gw_bayestar_url))
            return False
    else:
        msgs.warning('Failed to download: {}'.format(gw_bayestar_url))
        return False

    return True


def contours_from_gw_bayestar(file_name, credible_level=50.):
    r"""Given a bayestar.fits.gz HEALPix maps of a GW superevent it extracts contours at the given `credible_level`

    This is a wrapper around `ligo_skymap_contour`. The different contours are then counterclockwise oriented to be
    compatible with TAP and ASP queries.

    Args:
        file_name (`str`):
            Name of HEALPix maps
        credible_level (`float`):
            Probability level at which contours are returned

    Returns:
        contours (`list`):
           [RA,Dec] list defining the contours location. These are counterclowise oriented as seen on the sky from
           inside the sphere. The length of contours represent the number of not connected regions identified.
    """

    # Some checks
    assert isinstance(file_name, (str, np.str)), '{} is not a valid string'.format(file_name)
    if not checks.fits_file_is_valid(file_name):
        msgs.error('{} not a valid fits file'.format(file_name))
    assert isinstance(credible_level, (int, float, np.int, np.float)), '`credible_level` is not a float or an int'

    # Create temporary file where to store output from ligo_skymap_contour
    contour_tmp_file = '.'+file_name+'.tmp.json'
    if path.isfile(contour_tmp_file):
        remove(contour_tmp_file)

    # Use ligo_skymap_contour to compute contours
    ligo_skymap_contour(args=[file_name, '--output', contour_tmp_file, '--contour', str(credible_level)])

    # Parse the content of the resulting json file into a dict
    with open(contour_tmp_file, 'r') as jfile:
        jdata = jfile.read()
    contours_dict = json.loads(jdata)
    jfile.close()
    # cleaning up
    remove(contour_tmp_file)

    # Parse the resulting json to obtain a list of coordinates defining the vertices of the contours of each peak
    # contours is a list, each element of which contains the vertices of one contour encircling the desired significance
    contours = contours_dict['features'][0]['geometry']['coordinates']
    # Make sure that the orientation of the polygons on the sky is the correct one for the TAP queries,
    # i.e. counterclowise as seen on the sky from inside the sphere.
    for iii, contour in enumerate(contours):
        contours[iii] = _ensure_orientation(contour)

    # Quick summary:
    if len(contours)>0:
        msgs.info('Extracted the contours for {} regions at {} credible level'.format(len(contours), credible_level))
    else:
        msgs.info('No contours extracted at {} credible level'.format(credible_level))
        contours = None

    return contours


def _ensure_orientation(contour):
    r"""Set orientation of a polyong on the sky to be counterclockwise

    Take the input polygon and turn it into counterclockwise orientation as seen on the sky from inside the
    unitary sphere. This is done  through a `spherical_geometry.polygon.SphericalPolygon` object, which
    upon creation enforces the correct orientation independently of the one of the input (they call it clockwise,
    because it's see from outside).

    Args:
        contour (`list`):
            list of the vertices of a polygon

    Returns:
        contour_oriented (`list`)
            oriented list of the vertices of a polygon
    """

    # create the spherical polygon, which by construction has the right orientation on the sky
    contour_array = np.asarray(contour, dtype=np.float32)
    spherical_polygon = polygons.SphericalPolygon.from_radec(contour_array[:,0], contour_array[:,1], degrees=True)

    # reformat the correctly-oriented polygon to the same format as the input
    contour_oriented = []
    for radec in spherical_polygon.to_radec():
        for ra, dec in zip(radec[0], radec[1]):
            contour_oriented.append([ra, dec])

    return contour_oriented


def _array_split(array_in, thr):
    r"""Split an input array if two consecutive elements in x differ
    more than thr. The purpose is to deal with RAs folding at 360 degrees.
    """
    xx_in, yy_in = array_in[:,0], array_in[:,1]

    #Rearrange the input arrays to so that its first element is the most eastern one.
    ii = np.argmax(xx_in)
    xx_in = np.append(xx_in[ii:], xx_in[0:ii])
    yy_in = np.append(yy_in[ii:], yy_in[0:ii])

    dxx = np.diff(xx_in)
    split_indices = np.where(np.abs(dxx)>thr)
    return np.split(xx_in, split_indices[0]+1), np.split(yy_in, split_indices[0]+1)


def show_contours_from_gw_bayestar(file_name, contours=None, cmap='cylon', contours_color='C3',
                                   show_figure=True, save_figure=None):
    r"""Show sky credibility from the input healpix map and plot the contours created with `contours_from_gw_bayestar`

    Args:
        file_name:
        contours:
        cmap:
        contours_color:
        show_figure (`bool`):
        save_figure (`str` or `None`):

    Returns:

    """
    # Some checks
    assert isinstance(file_name, (str, np.str)), '{} is not a valid string'.format(file_name)
    if not checks.fits_file_is_valid(file_name):
        msgs.error('{} not a valid fits file'.format(file_name))

    if contours is None:
        msgs.warning("No contours defined, showing `credible_level=50` contours")
        plot_contours = contours_from_gw_bayestar(file_name, credible_level=50.)
    else:
        plot_contours = contours

    # ToDo find a more elegant solution for this.
    # ligo.skymap.tool.ligo_skymap_contour introduces problems with the matplotlib backend. This i a workaround to
    # get it running.
    importlib.reload(matplotlib)
    from matplotlib import pyplot as plt
    try:
        # check if the script runs on a notebook.
        # https://stackoverflow.com/questions/23883394/detect-if-python-script-is-run-from-an-ipython-shell-or-run-from
        # -the-command-li
        __IPYTHON__
    except NameError:
        # Try to get a working gui
        # https://stackoverflow.com/questions/39026653/programmatically-choose-correct-backend-for-matplotlib-on-mac-os-x
        gui_env = matplotlib.rcsetup.all_backends
        for gui in gui_env:
            try:
                matplotlib.use(gui, force=True)
                print(gui)
                break
            except:
                continue
    else:
        matplotlib.use('nbAgg')

    # Read map and get object name
    map_data, map_header = healpy.read_map(file_name, h=True, verbose=False, dtype=None)
    object_name_list = [value for name, value in map_header if name == 'OBJECT']
    if len(object_name_list)>0:
        object_name = object_name_list[0]
    else:
        object_name = 'GW event - {}'.format(file_name)

    # start the plot

    plt.figure(figsize=(10., 7.))
    ax = plt.axes([0.1, 0.1, 0.9, 0.9], projection='astro degrees mollweide')

    ax.grid(True)
    ax.imshow_hpx(map_data, cmap=cmap, visible=True, zorder=1)
    ax.text(0.15, 0.95, object_name, horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, bbox=dict(facecolor='orange', alpha=0.8), fontsize=24)
    # Convert the contour vertex coordinates from world to pixels and draw the contours
    w = wcs.WCS(ax.header)

    # Split a contour if two neighboring pixels are more distant than 10 degrees in RA
    split_step = (10 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / 360.)
    for contour_world in plot_contours:
        contour_pix = w.wcs_world2pix(contour_world, 0)
        x_contours, y_contours = _array_split(contour_pix, split_step)
        for x_contour, y_contour in zip(x_contours, y_contours):
            ax.plot(x_contour, y_contour, linewidth=2.5, color=contours_color, zorder=5)  # , marker='o')

    if save_figure is not None:
        plt.savefig(save_figure, dpi=200., format='pdf', bbox_inches='tight')

    if show_figure:
        plt.show()

    plt.close()

    # Bringing back the previously used matplotlib backend
    matplotlib.use(STARTING_MATPLOTLIB_BACKEND)





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
