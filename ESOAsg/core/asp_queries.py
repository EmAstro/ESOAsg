r"""Module to create and run ASP queries

The ESO Archive Science Portal ('ASP <http://archive.eso.org/scienceportal/home>'_) allows browsing and exploration
of archive content using an intuitive interactive user interface that supports iterative queries

"""

import webbrowser

from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.ancillary import checks

CONNECTORS = ['&']


def run_query(query, show_link=True, open_link=True):
    r"""Run the ASP query

    Args:
        query (str): url of the ASP query
        open_link (bool): open a link to the ASP page
        show_link (bool): show the link on the terminal

    Returns:
        None

    """
    if show_link:
        msgs.info('The ASP link is:\n {}\n'.format(query))
    if open_link:
        webbrowser.open(query)
    return


def base_url():
    r"""Return the url for ASP

    Returns:
        str: url

    """
    return str(default.get_value('eso_asp_url'))


def sort_by(sort_type='-obs_date', connector=None):
    r"""

    Args:
        sort_type (str): value for which sort the ASP query
        connector (str, optional): connector to be put in front of the the condition (e.g., '&')

    Returns:
        str: string that defines the sort type in a ASP query

    """
    if sort_type is None:
        return ''
    else:
        return '{}sort={}'.format(_get_connector(connector), sort_type)


def _get_connector(connector):
    r"""Check that the connector is valid

    Args:
        connector (str, optional): connector to be put in front of the the condition (e.g., '&')

    Returns:
        str: valid connector

    """
    if connector is None:
        clean_connector = ''
    elif connector in CONNECTORS:
        clean_connector = connector
    else:
        msgs.error('Not a valid connector for the ASP query. Possible values are: {}'.format(CONNECTORS))
        clean_connector = None
    return clean_connector


def condition_position(ra, dec, radius=None, connector=None):
    r"""Return condition for position

    Args:
        ra (any): string containing the RA of the source in degree
        dec (any): string containing the Dec of the source in degree
        radius (any): string containing the search radius in arcseconds
        connector (str, optional): connector to be put in front of the the condition (e.g., '&')

    Returns:
        str: string containing the `pos=` condition for ASP

    """
    if ra is None or dec is None:
        return ''
    condition_for_position = '{0}pos={1},{2}{3}'.format(_get_connector(connector),
                                                        checks.from_number_to_string(ra),
                                                        checks.from_number_to_string(dec),
                                                        condition_radius(radius))
    return condition_for_position


def condition_radius(radius, connector=None):
    r"""Return condition for radius

    Args:
        radius (any): string containing the search radius in arcseconds
        connector (str, optional): connector to be put in front of the the condition (e.g., '&')

    Returns:
        str: string containing the `pos=` condition for ASP

    """
    if radius is None:
        return ''
    condition_for_radius = '{0}pos={1}'.format(_get_connector(connector),
                                               checks.from_number_to_string(radius))
    return condition_for_radius


def from_polygons(polygons=None, open_link=False, show_link=False):
    if polygons is not None:
        for iii, polygon in enumerate(polygons):
            url = 'http://archive.eso.org/scienceportal/home?' + 'poly=' + polygon + '&sort=-obs_date'
            if show_link:
                msgs.info('ASP link to region N.{} is:\n {}\n'.format(str(iii + 1), url))
            if open_link:
                webbrowser.open(url)
