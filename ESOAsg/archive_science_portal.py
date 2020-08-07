from astropy import coordinates
from astropy.coordinates import ICRS

# from ESOAsg import msgs
from ESOAsg.ancillary import cleaning_lists
from ESOAsg.core import asp_queries


def query_from_radec(positions, radius=None, instruments=None, data_types=None, open_link=False, show_link=False):
    r"""Query the ESO `ASP service <http://archive.eso.org/scienceportal/home>`_ given a position

     The `positions` value (or list) needs to be given as an `astropy.coordinates.SkyCoord` object. For further detail
     see here: `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_

    Args:
        positions (astropy.coordinates.SkyCoord): coordinates (or list of coordinates) of the sky you want to query
        radius (float): search radius in arcseconds
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        open_link (bool): open a link to the ASP page
        show_link (bool): show the link on the terminal

    Returns:
        None

    """
    # Check input position
    positions_list = cleaning_lists.from_element_to_list(positions, element_type=coordinates.SkyCoord)
    # Working on instruments
    instruments_list = cleaning_lists.from_element_to_list(instruments, element_type=str)
    # Working on data_types
    data_types_list = cleaning_lists.from_element_to_list(data_types, element_type=str)
    # Run one query per position
    for iii, position in enumerate(positions_list):
        position.transform_to(ICRS)
        url = '{0}{1}{2}{3}{4}'.format(asp_queries.base_url(),
                                       asp_queries.condition_position(position.ra.degree, position.dec.degree, radius,
                                                                      connector=None),
                                       asp_queries.condition_instruments(instruments_list, connector='&'),
                                       asp_queries.condition_data_types(data_types_list, connector='&'),
                                       asp_queries.sort_by('-obs_date', connector='&'))
        asp_queries.run_query(url, show_link=show_link, open_link=open_link)

    return


def query_from_polygons(polygons=None, instruments=None, data_types=None, open_link=False, show_link=False):
    r"""Query the ESO `ASP service <http://archive.eso.org/scienceportal/home>`_ given a polygon

    The `polygons` value (or list) needs to be given as a string defining the location in the sky of the polygon
    with RA, Dec, separated by commas and with the first RA, Dec pair that matches the last one (to close the
    polygon)

    Args:
        polygons (list): list of `str` (or single `str`) containing the coordinates of the polygon in the sky you want
            to query
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        open_link (bool): open a link to the ASP page
        show_link (bool): show the link on the terminal

    Returns:
        None

    """
    # Check input position
    polygons_list = cleaning_lists.from_element_to_list(polygons, element_type=str)
    # Working on instruments
    instruments_list = cleaning_lists.from_element_to_list(instruments, element_type=str)
    # Working on data_types
    data_types_list = cleaning_lists.from_element_to_list(data_types, element_type=str)
    # Run one query per position
    for iii, polygon in enumerate(polygons_list):
        url = '{0}{1}{2}{3}{4}'.format(asp_queries.base_url(),
                                       asp_queries.condition_polygon(polygon, connector=None),
                                       asp_queries.condition_instruments(instruments_list, connector='&'),
                                       asp_queries.condition_data_types(data_types_list, connector='&'),
                                       asp_queries.sort_by('-obs_date', connector='&'))
        asp_queries.run_query(url, show_link=show_link, open_link=open_link)

    return
