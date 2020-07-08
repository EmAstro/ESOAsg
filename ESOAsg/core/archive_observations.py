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


def _run_query(query, verbose=False, remove_bytes=True, maxrec=default.get_value('maxrec')):
    r"""Run tap query and return result as a table

    Args:
        query (`str`):
            Query to be run
        verbose (`bool`):
            if set to `True` additional info will be displayed
        remove_bytes ('bool')
            if set to True, it converts all bytes entries to standard strings
        maxrec (`int`):
            Define the maximum number of entries that a single query can return

    Returns:
        result_from_query (`astropy.Table`):
            Result from the query to the TAP service
    """
    # Load tap service
    tapobs = _define_tap_service(verbose=False)
    if verbose:
        msgs.info('The query is: \n {} \n'.format(str(query)))
    # Obtaining query results and convert it to an astropy table
    result_from_query = tapobs.search(query=query, maxrec=maxrec).to_table()
    # removing bytes code from some columns:
    if remove_bytes:
        for column_name in result_from_query.colnames:
            result_from_query[column_name].data.data[:] = checks.from_bytes_to_string(result_from_query[
                                                                                          column_name].data.data)
    return result_from_query


def _query_obscore_base():
    r"""Create the base string for a query to `ivoa.ObsCore`

    Returns:
        query (`str`):
            Base for the ivoa.ObsCore query:
            SELECT
                target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max,
                dataproduct_type, instrument_name, obstech, abmaglim,
                proposal_id, obs_collection
            FROM
                ivoa.ObsCore
    """
    query_base = '''
            SELECT
                target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max, 
                dataproduct_type, instrument_name, obstech, abmaglim,
                proposal_id, obs_collection
            FROM
                ivoa.ObsCore'''
    return query_base


def _query_obscore_intersect_ra_dec(ra, dec, radius=None):
    r"""Create the where condition string for a query to `ivoa.ObsCore`

    Args:
        ra (`float`):
            RA of the target in degrees and in the ICRS system
        dec (`float`):
            Dec of the target in degrees and in the ICRS system
        radius (`float`):
            Search radius in arcsec. If set to `None` no radius will be considered in the INTERSECT condition
    Returns:
        query_intersect_ra_dec (`str`):
            String containing the WHERE INTERSECT condition for a query
    """
    if radius is None:
        query_intersect_ra_dec = '''
            WHERE
                INTERSECTS(POINT('ICRS',{},{}), s_region)=1'''.format(str(ra), str(dec))
    else:
        query_intersect_ra_dec = '''
            WHERE
                INTERSECTS(s_region,CIRCLE('ICRS',{},{},{}/3600.))=1'''.format(str(ra), str(dec), str(radius))
    return query_intersect_ra_dec


def _query_obscore_select_instruments(instruments_list):
    r"""Create condition string to select only specific instruments in `ivoa.ObsCore`

    Args:
        instruments_list (`list'):
            Limit the search to the selected list of instruments (e.g., `XSHOOTER`)
    Returns:
        query_select_instruments (`str`):
            String containing the `instrument_name=` condition for a query
    """
    if len(instruments_list) == 1:
        query_select_instruments = '''
            AND
                instrument_name='{}' '''.format(instruments_list[0])
    else:
        query_select_instruments = '''
            AND
                ('''
        for instrument_name in instruments_list:
            query_select_instruments = query_select_instruments + '''instrument_name='{}' OR '''.format(
                instrument_name)
        query_select_instruments = query_select_instruments[0:-4] + ')'
    return query_select_instruments


def _query_obscore_select_data_types(data_types_list):
    r"""Create condition string to select only specific dataproduct types in `ivoa.ObsCore`

    Args:
        data_types_list (`list'):
            Limit the search to the selected list of dataproduct types (e.g., `spectrum`)
    Returns:
        query_select_data_types (`str`):
            String containing the `dataproduct_type=` condition for a query
    """
    if len(data_types_list) == 1:
        query_select_data_types = '''
            AND
                dataproduct_type='{}' '''.format(data_types_list[0])
    else:
        query_select_data_types = '''
            AND
                ('''
        for data_type in data_types_list:
            query_select_data_types = query_select_data_types + '''dataproduct_type='{}' OR '''.format(data_type)
        query_select_data_types = query_select_data_types[0:-4] + ')'
    return query_select_data_types


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
        msgs.work('Downloading file {}. This may take some time.'.format(file_name + '.fits'))
        urllib.request.urlretrieve(download_url, filename=file_name + '.fits')
        msgs.info('File {} downloaded.'.format(file_name + '.fits'))


def query_from_radec(positions, radius=None, instruments=None, data_types=None, verbose=False,
                     maxrec=default.get_value('maxrec')):
    r"""Query the ESO TAP service given a position in RA and Dec.

     The `positions` value (or list) needs to be given as an `astropy.coordinates.SkyCoord` object.
    
    Args:
        positions (`astropy.coordinates.SkyCoord`):
            Coordinates (or list of coordinates) of the sky you want to query in the format of an
            `astropy.coordinates.SkyCoord` object. For further detail see here:
            `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_
        radius (`float`):
            Search radius you want to query in arcseconds. Note that in case `None` is given, the query will be
            performed with the `INTERSECT(POINT('',RA,Dec), s_region)` clause instead of the
            `INTERSECT(s_region,CIRCLE('',RA,Dec,radius/3600.))` one. See here for further examples:
            `tap obs examples <http://archive.eso.org/tap_obs/examples>`_
        instruments (`list`):
            Limit the search to the selected list of instruments (e.g., `XSHOOTER`)
        data_types (`list`):
            Limit the search to the selected types of data (e.g., `spectrum`)
        verbose (`bool`):
            if set to `True` additional info will be displayed
        maxrec (`int`):
            Define the maximum number of file that a single query can return from the ESO archive. The default values
            is set in the `default.txt` file.

    Returns:
        results_from_query (`list`):
            Results from the query in a list with the same length of the input position. Currently it contains:
            target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max, em_min, dataproduct_type, instrument_name,
            abmaglim, proposal_id, obs_collection
    """
    # Check inputs:
    # Working on positions
    if isinstance(positions, list):
        positions_list = positions
    else:
        positions_list = [positions]
    for position in positions_list:
        assert isinstance(position, coordinates.SkyCoord), r'Input positions not a SkyCoord object'
    # Working on radius
    if radius is not None:
        if isinstance(radius, int):
            radius = float(radius)
        else:
            assert isinstance(radius, float), r'Input radius is not a number'
    # Working on instruments
    if instruments is not None:
        if isinstance(instruments, list):
            instruments_list = instruments
        else:
            instruments_list = [instruments]
        for instrument in instruments_list:
            assert isinstance(instrument, str), r'Input instrument: {} not valid'.format(instrument)
    # Working on data_types
    if data_types is not None:
        if isinstance(data_types, list):
            data_types_list = data_types
        else:
            data_types_list = [data_types]
        for data_type in data_types_list:
            assert isinstance(data_type, str), r'Input data type: {} not valid'.format(data_type)

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
        ra, dec = np.float_(position.ra.degree), np.float_(position.dec.degree)
        msgs.work('Running query {} to the ESO archive (out of {} total)'.format(idx + 1, len(positions_list)))

        # Define query
        # base query:
        query = _query_obscore_base()
        # selection of the location:
        query = query + _query_obscore_intersect_ra_dec(ra, dec, radius=radius)
        # selection of the instrument(s)
        if instruments is not None:
            query = query + _query_obscore_select_instruments(instruments_list)
        # selection of the data_type(s)
        if data_types is not None:
            query = query + _query_obscore_select_data_types(data_types_list)

        # running query and append results to the list
        result_from_query = _run_query(query, verbose=verbose, remove_bytes=True, maxrec=maxrec)

        if len(result_from_query) < 1:
            msgs.warning('No data has been retrieved')
        else:
            msgs.info('A total of {} entries has been retrieved'.format(len(result_from_query)))
            if verbose:
                msgs.info('For the following instrument:')
                for inst_name in np.unique(result_from_query['instrument_name'].data):
                    msgs.info(' - {}'.format(inst_name))

        results_from_query.append(result_from_query)
    return results_from_query


def query_ASP_from_radec(positions, radius=None, open_link=False, show_link=False):
    r"""Query the ESO ASP service given a position in RA and Dec.

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
        open_link (`bool`):
            open a link to the ASP page
        show_link (`bool`):
            show the link on the terminal

    """
    # Check inputs:
    # Working on positions
    if isinstance(positions, list):
        positions_list = positions
    else:
        positions_list = [positions]
    for position in positions_list:
        assert isinstance(position, coordinates.SkyCoord), r'Input positions not a SkyCoord object'
    # Working on radius
    if radius is not None:
        if isinstance(radius, int):
            radius = float(radius)
        else:
            assert isinstance(radius, float), r'Input radius is not a number'

    for iii, position in enumerate(positions_list):
        position.transform_to(ICRS)
        radec_string = np.str(np.float_(position.ra.degree)) + ',' + np.str(np.float_(position.dec.degree))
        radius_string = np.str(radius / 3600.)
        url = 'http://archive.eso.org/scienceportal/home?' + 'pos=' + radec_string + '&r=' + radius_string + \
              '&sort=-obs_date'
        if show_link:
            msgs.info('ASP link to region N.{} is:\n {}\n'.format(np.str(iii + 1), url))
        if open_link:
            webbrowser.open(url)


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
            list_of_outputs = [output + text_file for output in list_of_files]
    elif isinstance(text_file, list):
        if len(list_of_files) == len(text_file):
            list_of_outputs = text_file
        else:
            list_of_outputs = [files + '.hdr' for files in list_of_files]
    else:
        list_of_outputs = [files + '.hdr' for files in list_of_files]

    # Downloading headers
    for file_name, file_out in zip(list_of_files, list_of_outputs):
        if os.path.isfile(file_out):
            msgs.warning('Overwriting existing text file: {}'.format(file_out))
            os.remove(file_out)
        url_for_header = archive_url + 'hdr?DpId=' + file_name
        response_url = requests.get(url_for_header, allow_redirects=True)
        # Removing html from text
        header_txt = response_url.text.split('<pre>')[1].split('</pre>')[0]
        if not header_txt.startswith('No info found for'):
            file_header = open(file_out, 'w')
            for line in header_txt.splitlines():
                file_header.write(line + '\n')
            file_header.close()
            msgs.info('Header successfully saved in: {}'.format(file_out))
        else:
            msgs.warning('{} is not present in the ESO archive'.format(file_name))

    return


def query_ASP_from_polygons(polygons=None, open_link=False, show_link=False):
    if polygons is not None:
        for iii, polygon in enumerate(polygons):
            url = 'http://archive.eso.org/scienceportal/home?' + 'poly=' + polygon + '&sort=-obs_date'
            if show_link:
                msgs.info('ASP link to region N.{} is:\n {}\n'.format(np.str(iii + 1), url))
            if open_link:
                webbrowser.open(url)


def query_TAP_from_polygons(polygons=None, merge=False, instrument=None, maxrec=default.get_value('maxrec'),
                            verbose=False):
    tapobs = dal.tap.TAPService(default.get_value('eso_tap_obs'))
    msgs.info('Querying the ESO TAP service at:')
    msgs.info('{}'.format(str(default.get_value('eso_tap_obs'))))

    results_from_query = []

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
                results_from_query.append(result_from_query)

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
            results_from_query.append(result_from_query)
        return results_from_query

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
        polygon = ' '.join(['%.4f, %.4f,' % (ra, dec) for ra, dec in contour_clean])[
                  :-1]  # Remove the last character, which is an unwanted extra comma
        polygons.append(polygon)

    return polygons
