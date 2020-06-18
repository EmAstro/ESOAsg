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
from astropy import coordinates
from astropy.coordinates import ICRS
import requests
import webbrowser


from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.ancillary import checks


def _define_tap_service(verbose=False):
    r"""Load tap service from defaults

    The TAP service for raw. reduced, and ambient data is defined in `ESOAsg\default.txt` as `eso_tap_obs`

    Args:
        verbose (`bool`):
            if set to `True` additional info will be displayed

    Returns:
        tapcat (`pyvo.dal.tap.TAPService`)
            TAP service that will be used for the queries
    """
    if verbose:
        msgs.info('Querying the ESO TAP service at:')
        msgs.info('{}'.format(str(default.get_value('eso_tap_obs'))))
    tapobs = dal.tap.TAPService(default.get_value('eso_tap_obs'))
    return tapobs


def _run_query(query, verbose=False, maxrec=default.get_value('maxrec')):
    r"""Run tap query and return result as a table

    Args:
        query (`str`):
            Query to be run
        verbose (`bool`):
            if set to `True` additional info will be displayed
        maxrec (`int`, `None`):
            Define the maximum number of entries that a single query can return

    Returns:
        result_from_query (`astropy.Table`):
            Result from the query to the TAP service
    """
    # Load tap service
    tapobs = _define_tap_service(verbose=False)
    if verbose:
        msgs.info('The query is:')
        msgs.info('{}'.format(str(query)))
    # Obtaining query results and convert it to an astropy table
    result_from_query = tapobs.search(query=query, maxrec=maxrec).to_table()
    return result_from_query


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


def query_from_radec(positions, radius=None, instruments=None, verbose=False, maxrec=default.get_value('maxrec')):
    r"""Query the ESO TAP service given a position in RA and Dec.

     The `positions` value (or list) needs to be given as an `astropy.coordinates.SkyCoord` object.
    
    Args:
        positions (`astropy.coordinates.SkyCoord`):
            Coordinates (or list of coordinates) of the sky you want to query in the format of an
            `astropy.coordinates.SkyCoord` object. For further detail see here:
            `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_
        radius (`float`):
            Search radius you want to query in arcseconds. Note that in case `None` is given, the query will be
            performed with the `CONTAINS(POINT('',RA,Dec), s_region)` clause instead of the
            `CONTAINS(s_region,CIRCLE('',RA,Dec,radius/3600.))` one. See here for further examples:
            `tap obs examples <http://archive.eso.org/tap_obs/examples>`_
        instruments (`list'):
            Limit the search to the selected list of instruments
        verbose (`bool`):
            if set to `True` additional info will be displayed
        maxrec (`int`):
            Define the maximum number of file that a single query can return from the ESO archive. You probably never
            need this. By default is set by the `default.txt` file.

    Returns:
        results_from_query (`pyvo.dal.tap.TAPResults`):
            Result from the query. Currently it contains: target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max,
            em_min, dataproduct_type, instrument_name, abmaglim, proposal_id

    """
    # Check inputs:
    # Working on positions
    if isinstance(positions, list):
        positions_list = positions
    else:
        positions_list = [positions]
    for position in positions_list:
        assert isinstance(position, coordinates.SkyCoord), r'Input positions not a SkyCoord object'
    # Working on instruments
    if instruments is not None:
        if isinstance(instruments, list):
            instruments_list = instruments
        else:
            instruments_list = [instruments]
        for instrument in instruments_list:
            assert isinstance(instrument, str), r'Input instrument: {} not valid'.format(instrument)
    # Working on radius
    if radius is not None:
        if isinstance(radius, int):
            radius = float(radius)
        else:
            assert isinstance(radius, float), r'Input radius is not a number'

    # Running over all positions
    if verbose:
        how_many_positions = len(positions_list)
        if how_many_positions > 1:
            msgs.work('Exploring ESO archive around {} locations in the sky'.format(how_many_positions))
        else:
            msgs.work('Exploring ESO archive around the input location in the sky')

    results_from_query = []

    for position, idx in zip(positions_list, range(len(positions_list))):
        position.transform_to(ICRS)
        ra, dec = np.float32(position.ra.degree), np.float32(position.dec.degree)
        msgs.working('Running query {} out of {} to the ESO'.format(idx+1,len(positions_list)))
        # Define query
        if radius is None:
            query = '''
                    SELECT
                        target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max, 
                        dataproduct_type, instrument_name, abmaglim, proposal_id
                    FROM
                        ivoa.ObsCore
                    WHERE
                        CONTAINS(POINT('ICRS',{},{}), s_region)=1
                    '''.format(ra, dec)
        else:
            query = '''
                    SELECT
                        target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max, 
                        dataproduct_type, instrument_name, abmaglim, proposal_id
                    FROM
                        ivoa.ObsCore
                    WHERE
                        CONTAINS(s_region,CIRCLE('ICRS',{},{},{}/3600.))=1
                    '''.format(ra, dec, radius)
        if instruments is not None:
            if len(instruments_list) == 1:
                instruments_selection = '''
                                        AND
                                            instrument_name="{}"
                                        '''.format(instruments_list[0])
            else:
                instruments_selection = '''
                                        AND
                                            ('''
                for instrument_name in instruments_list:
                    instruments_selection = instruments_selection + '''instrument_name="{}" OR '''.format(
                        instrument_name)
                instruments_selection = instruments_selection[0:-3] + ')'
            query = '\n'.join([query, instruments_selection])
        msgs.info('The query is:')
        msgs.info('{}'.format(str(query)))
        results_from_query.append(_run_query(query, verbose=verbose))
    return results_from_query

    """
    
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
    """

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


def query_ASP_from_polygons(polygons=None, open_link=False, show_link=False):
    if polygons is not None:
        for iii, polygon in enumerate(polygons):
            url = 'http://archive.eso.org/scienceportal/home?' + 'poly='+polygon + '&sort=-obs_date'
            if show_link:
                msgs.info('ASP link to region N.{} is:\n {}\n'.format(np.str(iii+1),url))
            if open_link:
                webbrowser.open(url)


def query_TAP_from_polygons(polygons=None, merge=False, instrument=None, maxrec=default.get_value('maxrec'),
                            verbose=False):

    tapobs = dal.tap.TAPService(default.get_value('eso_tap_obs'))
    msgs.info('Querying the ESO TAP service at:')
    msgs.info('{}'.format(str(default.get_value('eso_tap_obs'))))

    result = []

    if polygons is not None:
        polygon_union = ''
        for iii, polygon in enumerate(polygons):
            polygon_union += """intersects(s_region, POLYGON('', """ + polygon + """)) = 1 OR """
            if not merge:
                query = """SELECT
                               target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max, 
                               dataproduct_type, instrument_name, abmaglim, proposal_id
                           FROM
                               ivoa.ObsCore
                           WHERE
                               intersects(s_region, POLYGON('', """ + polygon + """)) = 1 """
                if instrument is not None:
                    instrument_selection = str(
                        """                     AND instrument_name='{}'""".format(str(instrument)))
                    query = '\n'.join([query, instrument_selection])
                if verbose:
                    msgs.info('The query is:')
                    msgs.info('{}'.format(str(query)))

                # Obtaining query results
                result_from_query = tapobs.search(query=query, maxrec=maxrec)
                if len(result_from_query) < 1:
                    msgs.warning('No data has been retrieved')
                else:
                    msgs.info('A total of {} entries has been retrieved for polygon N.{}'.format(len(
                        result_from_query), iii))
                    msgs.info('For the following instrument:')
                    for inst_name in np.unique(result_from_query['instrument_name'].data):
                        msgs.info(' - {}'.format(inst_name.decode("utf-8")))
                    if verbose:
                        result_from_query.to_table().pprint(max_width=-1)
                result.append(result_from_query)

        polygon_union = polygon_union[:-4]
        if merge:
            query = """SELECT
                           target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max, 
                           dataproduct_type, instrument_name, abmaglim, proposal_id
                       FROM
                           ivoa.ObsCore
                       WHERE
                            (""" + polygon_union + """)"""
            if instrument is not None:
                instrument_selection = str(
                    """                     AND instrument_name='{}'""".format(str(instrument)))
                query = '\n'.join([query, instrument_selection])
            if verbose:
                msgs.info('The query is:')
                msgs.info('{}'.format(str(query)))
            # Obtaining query results
            result_from_query = tapobs.search(query=query, maxrec=maxrec)
            if len(result_from_query) < 1:
                msgs.warning('No data has been retrieved')
            else:
                msgs.info('A total of {} entries has been retrieved for polygon N.{}'.format(len(
                    result_from_query), iii))
                msgs.info('For the following instrument:')
                for inst_name in np.unique(result_from_query['instrument_name'].data):
                    msgs.info(' - {}'.format(inst_name.decode("utf-8")))
                if verbose:
                    result_from_query.to_table().pprint(max_width=-1)
            result.append(result_from_query)
        return result

        '''
        query="""SELECT count(*) from ivoa.ObsCore
                 WHERE (""" + polygon_union + """)"""
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
        '''


def contours_to_polygons(contours, max_vertices=30):
    r"""

    Args:
        contours:
        max_vertices (`int`):

    Returns:
        polygons:

    """
    polygons = []
    for iii, contour in enumerate(contours):
        if len(contour) > max_vertices:
            contour_clean = contour[0: len(contour): int(len(contour) / max_vertices + 1)]
        else:
            contour_clean = contour
        # Construct the polygon as a string in the right format
        polygon = ' '.join(['%.4f, %.4f,' % (ra, dec) for ra, dec in contour_clean])[:-1] # Remove the last character, which is an unwanted extra comma
        polygons.append(polygon)

    return polygons
