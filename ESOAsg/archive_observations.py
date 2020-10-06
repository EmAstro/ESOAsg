from astropy import coordinates
from astropy.coordinates import ICRS

import numpy as np
import urllib

from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.core import tap_queries
from ESOAsg.queries import query_observations
from ESOAsg.ancillary import checks
from ESOAsg.ancillary import cleaning_lists


def query_from_radec(positions=None, radius=None, instruments=None, data_types=None, verbose=False,
                     maxrec=None):
    r"""Query the ESO archive for data at a given position in RA and Dec

    The `positions` value (or list) needs to be given as an
    `astropy.coordinates.SkyCoord <https://docs.astropy.org/en/stable/coordinates/>`_ object.

    The output is in an (list of) `astropy.table` with columns defined in: `core.tap_queries.COLUMNS_FROM_OBSCORE`

    .. note::
        In case you are querying radius=`None` is set, the query will performed with:
        `INTERSECT(POINT('',RA,Dec), s_region)`
        instead of:
        `INTERSECT(s_region,CIRCLE('',RA,Dec,radius/3600.))`.
        See here for further examples: `tap obs examples <http://archive.eso.org/tap_obs/examples>`_


    Args:
        positions (astropy.coordinates.SkyCoord): coordinates (or list of coordinates) of the sky you want to query
        radius (float, optional): search radius in arcseconds
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        verbose (bool): if set to `True` additional info will be displayed
        maxrec (int, optional): define the maximum number of entries that a single query can return. If it is `None` the
            value is set by the limit of the service.

    Returns:
        any: results from the queries

    """
    # Check inputs:
    # Working on positions
    positions_list = cleaning_lists.from_element_to_list(positions, element_type=coordinates.SkyCoord)
    # Working on radius
    if radius is not None:
        if isinstance(radius, int):
            radius = float(radius)
        else:
            assert isinstance(radius, float), r'Input radius is not a number'
    # Working on instruments
    instruments_list = cleaning_lists.from_element_to_list(instruments, element_type=str)
    # Working on data_types
    data_types_list = cleaning_lists.from_element_to_list(data_types, element_type=str)

    if verbose:
        how_many_positions = len(positions_list)
        if how_many_positions > 1:
            msgs.work('Exploring ESO archive around {} locations in the sky'.format(how_many_positions))
        else:
            msgs.work('Exploring ESO archive around the input location in the sky')

    # Running over all positions
    results_from_query = []
    for idx, position in enumerate(positions_list):
        position.transform_to(ICRS)
        ra, dec = np.float_(position.ra.degree), np.float_(position.dec.degree)
        msgs.work('Running query {} to the ESO archive (out of {} total)'.format(idx + 1, len(positions_list)))
        # Define query
        query = "{0}{1}{2}{3}".format(tap_queries.create_query_obscore_base(),
                                      tap_queries.condition_intersects_ra_dec(ra, dec, radius=radius),
                                      tap_queries.condition_instruments_like(instruments_list),
                                      tap_queries.condition_data_types_like(data_types_list))
        # instantiate ESOCatalogues
        query_for_observations = query_observations.ESOObservations(query=query, type_of_query='sync', maxrec=maxrec)
        # running query and append results to the list
        if verbose:
            query_for_observations.print_query()
        # Obtaining query results
        query_for_observations.run_query(to_string=True)
        result_from_query = query_for_observations.get_result_from_query()
        if len(result_from_query) < 1:
            msgs.warning('No data has been retrieved')
        else:
            msgs.info('A total of {} entries has been retrieved'.format(len(result_from_query)))
            if verbose:
                msgs.info('For the following instrument:')
                for inst_name in np.unique(result_from_query['instrument_name'].data):
                    msgs.info(' - {}'.format(inst_name))

        results_from_query.append(result_from_query)

    # Returning results
    if len(results_from_query) == 0:
        return None
    elif len(results_from_query) == 1:
        return results_from_query[0]
    else:
        return results_from_query


def query_from_polygons(polygons, instruments=None, data_types=None, verbose=False, maxrec=None):
    r"""Query the ESO archive for data at a area in the sky defined by a polygon

    The `polygons` value (or list) needs to be given as a string defining the location in the sky of the polygon
    with RA, Dec, separated by commas and with the first RA, Dec pair that matches the last one (to close the
    polygon)

    The output is in an (list of) `astropy.table` with columns defined in: `core.tap_queries.COLUMNS_FROM_OBSCORE`


    Args:
        polygons (list): ist of `str` (or single `str`) containing the coordinates of the polygon in the sky you want
            to query
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        verbose (bool): if set to `True` additional info will be displayed
        maxrec (int, optional): define the maximum number of entries that a single query can return. If it is `None` the
            value is set by the limit of the service.

    Returns:
        any: results from the queries

    """
    # Check inputs:
    # Working on polygons
    polygons_list = cleaning_lists.from_element_to_list(polygons, element_type=str)
    # Working on instruments
    instruments_list = cleaning_lists.from_element_to_list(instruments, element_type=str)
    # Working on data_types
    data_types_list = cleaning_lists.from_element_to_list(data_types, element_type=str)

    if verbose:
        how_many_polygons = len(polygons_list)
        if how_many_polygons > 1:
            msgs.work('Exploring ESO archive around {} locations in the sky'.format(how_many_polygons))
        else:
            msgs.work('Exploring ESO archive around the input location in the sky')

    # Running over all positions
    results_from_query = []
    for idx, polygon in enumerate(polygons_list):
        msgs.work('Running query {} to the ESO archive (out of {} total)'.format(idx + 1, len(polygons_list)))
        # Define query
        query = "{0}{1}{2}{3}".format(tap_queries.create_query_obscore_base(),
                                      tap_queries.condition_intersects_polygon(polygon),
                                      tap_queries.condition_instruments_like(instruments_list),
                                      tap_queries.condition_data_types_like(data_types_list))
        # instantiate ESOCatalogues
        query_for_observations = query_observations.ESOObservations(query=query, type_of_query='sync', maxrec=maxrec)
        # running query and append results to the list
        if verbose:
            query_for_observations.print_query()
        # Obtaining query results
        query_for_observations.run_query(to_string=True)
        result_from_query = query_for_observations.get_result_from_query()
        if len(result_from_query) < 1:
            msgs.warning('No data has been retrieved')
        else:
            msgs.info('A total of {} entries has been retrieved (with maxrec={})'.format(len(result_from_query),
                                                                                         maxrec))
            msgs.info('For the following instrument:')
            for inst_name in np.unique(result_from_query['instrument_name'].data):
               msgs.info(' - {}'.format(inst_name))

        results_from_query.append(result_from_query)

    # Returning results
    if len(results_from_query) == 0:
        return None
    elif len(results_from_query) == 1:
        return results_from_query[0]
    else:
        return results_from_query


def download(dp_ids, min_disk_space=float(default.get_value('min_disk_space'))):
    r"""Given a filename in the ADP format, the code download the file from the
    `ESO archive <http://archive.eso.org>`_

    Args:
        dp_ids (any): list data product ID (or single product ID) to be downloaded
        min_disk_space (float): the file will be downloaded only if there is this amount of space (in Gb) free on the
            disk

    Returns:
        None

    """
    # Check for disk space
    checks.check_disk_space(min_disk_space=min_disk_space)
    # Cleaning list
    dp_ids_list = cleaning_lists.from_element_to_list(cleaning_lists.from_bytes_to_string(dp_ids), element_type=str)
    for dp_id in dp_ids_list:
        # Given a dp_id of a public file, the link to download it is constructed as follows:
        download_url = 'http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{}&eso_download=file'.format(dp_id)
        msgs.work('Retrieving file {}.fits'.format(dp_id))
        urllib.request.urlretrieve(download_url, filename=dp_id + '.fits')
        msgs.info('File {}.fits downloaded'.format(dp_id))
