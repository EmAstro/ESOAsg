from astropy import coordinates
from astropy.coordinates import ICRS

# from ESOAsg import msgs
from ESOAsg.ancillary import checks
from ESOAsg.core import asp_queries


def query_from_radec(positions, radius=None, open_link=False, show_link=False):
    r"""Query the ESO `ASP service <http://archive.eso.org/scienceportal/home>`_ given a position

     The `positions` value (or list) needs to be given as an `astropy.coordinates.SkyCoord` object. For further detail
     see here: `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_

    Args:
        positions (astropy.coordinates.SkyCoord): coordinates (or list of coordinates) of the sky you want to query
        radius (float): search radius in arcseconds
        open_link (bool): open a link to the ASP page
        show_link (bool): show the link on the terminal

    Returns:
        None

    """
    # Check input position
    positions_list = checks.from_element_to_list(positions, element_type=coordinates.SkyCoord)
    # Run one query per position
    for iii, position in enumerate(positions_list):
        position.transform_to(ICRS)
        url = '{0}{1}{2}'.format(asp_queries.base_url(),
                                 asp_queries.condition_position(position.ra.degree, position.dec.degree, radius,
                                                                connector=None),
                                 asp_queries.sort_by('-obs_date', connector='&'))
        asp_queries.run_query(url, show_link=show_link, open_link=open_link)

    return
