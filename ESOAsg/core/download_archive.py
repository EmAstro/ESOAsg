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

# import os
# import sys
import urllib
import numpy as np
import os

from pyvo import dal
from astropy.coordinates import ICRS
import requests


from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.ancillary import checks


def download(dp_id, min_disk_space=np.float32(default.get_value('min_disk_space'))):
    r"""Given a filename in the ADP format, the code download the file from the
    `ESO archive <http://archive.eso.org>`_

    ..note::
        if dp_id is not a `numpy.str`, a WARNING message will be raised and the content of `dp_id` will be
        converted into a string.

    Args:
        dp_id (`numpy.str`):
            Data product ID to be downloaded.
        min_disk_space (`numpy.float`):
            The file will be downloaded only if there is this amount of space (in Gb) free on the disk.
            By default is set by the `default.txt` file.

    Returns:
        This downloads a fits ADP file with the same name of the input.
    """

    # Check for disk space
    checks.check_disk_space(min_disk_space=min_disk_space)

    for file_name in dp_id:
        # if the file name is in byte, this decode it.
        if not isinstance(file_name, str):
            msgs.warning('The content of dp_id is not in a string format.')
            msgs.warning('The code is trying to fix this.')
            if isinstance(file_name, bytes):
                file_name = np.str(file_name.decode("utf-8"))
                msgs.warning('Converted to {}.'.format(type(file_name)))
            else:
                msgs.error('Unable to understand the format of the dp_id entry: {}'.format(type(file_name)))

        # Given a dp_id of a public file, the link to download it is constructed as follows:
        download_url = 'http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{}&eso_download=file'.format(
            str(file_name))
        msgs.work('Downloading file {}. This may take some time.'.format(file_name+'.fits'))
        urllib.request.urlretrieve(download_url, filename=file_name + '.fits')
        msgs.info('File {} downloaded.'.format(file_name + '.fits'))


def query_from_radec(position, radius=None, instrument=None, maxrec=default.get_value('maxrec')):
    r"""
    Query to the ESO TAP service (link defined in `ESOAsg\default.txt`) for a specific location. The
    `position` needs to be given as an `astropy.coordinates.SkyCoord` object.
    
    Args:
        position (`astropy.coordinates.SkyCoord`):
            Coordinates of the sky you want to query in the format of an `astropy.coordinates.SkyCoord` object. Note
            that at the moment it works for one target at the time. For further detail see here:
            `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_
        radius (`numpy.float`):
            Search radius you want to query in arcseconds. Note that in case `None` is given, the query will be
            performed with the `CONTAINS(POINT('',RA,Dec), s_region)` clause instead of the
            `CONTAINS(s_region,CIRCLE('',RA,Dec,radius/3600.))` one. See here for further examples:
            `tap obs examples <http://archive.eso.org/tap_obs/examples>`_
        instrument (`str'):
            Limit the search to the selected instrument
            ToDo EMA:
            At the moment this works for a single instrument. Lists of instruments will be added in the
            future.
        maxrec (`numpy.int`):
            Define the maximum number of file that a single query can return from the ESO archive. You probably never
            need this. By default is set by the `default.txt` file.

    Returns:
        result_from_query (`pyvo.dal.tap.TAPResults`):
            Result from the query. Currently it contains: target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max,
            em_min, dataproduct_type, instrument_name, abmaglim, proposal_id

    """

    # Define TAP SERVICE
    tapobs = dal.tap.TAPService(default.get_value('eso_tap_obs'))
    msgs.info('Querying the ESO TAP service at:')
    msgs.info('{}'.format(str(default.get_value('eso_tap_obs'))))

    # ToDo EMA
    # This should be more flexible and take lists/arrays as input
    position.transform_to(ICRS)
    if not np.isscalar(position.ra.degree):
        if np.ndim(position.ra.degree) > 1:
            msgs.warning('The position should refer to a single pointing')
            msgs.warning('only the first location is taken into account')
            msgs.warning('multi pointing will be available in a future release')
        # Converting from array to number
        ra, dec = np.float32(position.ra.degree[0]), np.float32(position.dec.degree[0])
    else:
        ra, dec = np.float32(position.ra.degree), np.float32(position.dec.degree)

    # Define query
    if radius is not None:
        if not np.isscalar(radius):
            radius_scalar = np.float32(radius[0])
        else:
            radius_scalar = np.float32(radius)
        query = """SELECT
                      target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max, 
                      dataproduct_type, instrument_name, abmaglim, proposal_id
                   FROM
                     ivoa.ObsCore
                   WHERE
                     CONTAINS(s_region,CIRCLE('ICRS',{},{},{}/3600.))=1""".format(ra, dec, radius_scalar)
    else:
        query = """SELECT
                     target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max, 
                     dataproduct_type, instrument_name, abmaglim, proposal_id
                   FROM
                     ivoa.ObsCore
                   WHERE
                     CONTAINS(POINT('ICRS',{},{}), s_region)=1""".format(ra, dec)
    if instrument is not None:
        instrument_selection = str("""                     AND instrument_name='{}'""".format(str(instrument)))
        query = '\n'.join([query, instrument_selection])
    msgs.info('The query is:')
    msgs.info('{}'.format(str(query)))

    # Obtaining query results
    result_from_query = tapobs.search(query=query, maxrec=maxrec)
    if len(result_from_query) < 1:
        msgs.warning('No data has been retrieved')
    else:
        msgs.info('A total of {} entries has been retrieved'.format(len(result_from_query)))
        msgs.info('For the following instrument:')
        for inst_name in np.unique(result_from_query['instrument_name'].data):
            msgs.info(' - {}'.format(inst_name.decode("utf-8")))
    return result_from_query


def get_header_from_archive(file_id, text_file=None):  # written by Ema. 04.03.2020
    r"""Given a file ID the macro download the corresponding header.

    Args:
        file_id (`str`):
            ESO file ID for which the header will be downloaded
        text_file (`str`):
            text file where the header will be downloaded. If `None` it will it will be set to the same
            string `file_id` but with a `.hdr` extension.

    """

    # checks for connection to ESO archive
    archive_url = default.get_value('eso_archive_url')
    if not checks.connection_to_website(archive_url, timeout=1):
        msgs.error('Cannot connect to the ESO archive website:\n {}'.format(archive_url))

    # checks for file id
    assert isinstance(file_id, list) or isinstance(file_id, (str, np.str)), 'file_id needs to be a str or a list'
    if isinstance(file_id, str):
        list_of_files = [file_id]
    else:
        list_of_files = file_id
    list_of_files = [files if not files.endswith('.fits') else files.replace('.fits', '') for files in list_of_files]

    # checks for text_file
    assert isinstance(text_file, list) or isinstance(text_file, (str, np.str)) or \
           isinstance(text_file, (type(None), bytes)), 'text_file needs to be a str or a list'
    if isinstance(text_file, str):
        if len(list_of_files) == 1:
            list_of_outputs = [text_file]
        else:
            list_of_outputs = [output+text_file for output in list_of_files]
    elif isinstance(text_file, list):
        if len(list_of_files) == len(text_file):
            list_of_outputs = text_file
        else:
            list_of_outputs = [files + '.hdr' for files in list_of_files]
    else:
        list_of_outputs = [files+'.hdr' for files in list_of_files]

    # Downloading headers
    for file_name, file_out in zip(list_of_files, list_of_outputs):
        if os.path.isfile(file_out):
            msgs.warning('Overwriting existing text file: {}'.format(file_out))
            os.remove(file_out)
        url_for_header = archive_url+'hdr?DpId='+file_name
        response_url = requests.get(url_for_header, allow_redirects=True)
        # Removing html from text
        header_txt = response_url.text.split('<pre>')[1].split('</pre>')[0]
        if not header_txt.startswith('No info found for'):
            file_header = open(file_out, 'w')
            for line in header_txt.splitlines():
                file_header.write(line+'\n')
            file_header.close()
            msgs.info('Header successfully saved in: {}'.format(file_out))
        else:
            msgs.warning('{} is not present in the ESO archive'.format(file_name))

    return
