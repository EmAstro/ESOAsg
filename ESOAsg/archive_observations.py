from astropy import coordinates

import numpy as np

from ESOAsg import msgs
from ESOAsg.ancillary import checks


def query_from_radec(positions, radius=None, instruments=None, data_types=None, verbose=False,
                     maxrec=None):
    r"""Query the ESO archive for data at a given a position in RA and Dec.

     The `positions` value (or list) needs to be given as an
     `astropy.coordinates.SkyCoord <https://docs.astropy.org/en/stable/coordinates/>`_ object.

    .. note::
        In case radius=`None` is set, the query will performed with: `INTERSECT(POINT('',RA,Dec), s_region)` instead of:
        `INTERSECT(s_region,CIRCLE('',RA,Dec,radius/3600.))`. See here for further examples: `tap obs examples
        <http://archive.eso.org/tap_obs/examples>`_


    Args:
        positions (astropy.coordinates.SkyCoord): coordinates (or list of coordinates) of the sky you want to query
        radius (float, optional): Search radius in arcseconds
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        verbose (bool): if set to `True` additional info will be displayed
        maxrec (int, optional): define the maximum number of entries that a single query can return. If it is `None` the
            value is set by the limit of the service.

    Returns:
        list: Results from the query in a list with the same length of the input position. Currently it contains:
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
    instruments_list = checks.from_string_to_list(instruments)

    # Working on data_types
    data_types_list = checks.from_string_to_list(data_types)

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