r"""
get_data_from_radec
===================
Script to collect data from the ESO archive given RA and Dec

.. topic:: Inputs:

    - **ra_degree** - RA of the target in degree [J2000]
    - **dec_degree** - Dec of the target in degree [J2000]
    - **radius** - Search cone radius in arcsec
    - **instruments** - Limit the search to the selected ESO instruments
    - **data_types** - Limit the search to the selected data type
    - **maxrec** - Maximum number of files you would like to retrieve
"""

import argparse
import sys
from ESOAsg import __version__
from ESOAsg import default

EXAMPLES = str(r"""EXAMPLES:""" + """\n""" + """\n""" +
               r"""Get images that intersect a circle of 5 arcseconds around the position: """ +
               r"""RA=15.054250 and Dec=-28.048833 """ + """\n""" +
               r"""`_changed` as ending (and update the `CHECKSUM`): """ + """\n""" +
               r""">>> get_data_from_radec """ +
               r"""-ra 15.054250 """ +
               r"""-dec -28.048833 """ +
               r"""-r 5. """ +
               r"""-dt image """ + """\n""" +
               r""" """)


def parser(options=None):
    parser = argparse.ArgumentParser(
        description=r"""Collect data from the ESO archive given RA and Dec in degrees.""" + """\n""" + """\n""" +
                    r"""This uses ESOAsg version {:s}""".format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument("-ra", "--ra_degree", nargs="+", type=float, default=None,
                        help=r"RA of the target in degree [J2000]")
    parser.add_argument("-dec", "--dec_degree", nargs="+", type=float, default=None,
                        help=r"Dec of the target in degree [J2000]")
    parser.add_argument("-r", "--radius", nargs="+", type=float, default=None,
                        help=r"Search cone radius in arcsec")
    parser.add_argument("-i", "--instruments", nargs="+", type=str, default=None,
                        help="Limit the search to the selected ESO instruments")
    parser.add_argument("-dt", "--data_types", nargs="+", type=str, default=None,
                        help="Limit the search to the selected data type")
    parser.add_argument("-m", "--maxrec", nargs="+", type=int, default=default.get_value('maxrec'),
                        help="Maximum number of files you would like to retrieve")
    parser.add_argument("-v", "--version", action="version", version=__version__)

    # allows for negative values in dec
    for i in sys.argv[1:]:
        if i.startswith('--dec_') or i.startswith('-dec'):
            sys.argv[sys.argv.index(i) + 1] = str(float(sys.argv[sys.argv.index(i) + 1]))

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    from astropy import coordinates
    from astropy import units as u

    from ESOAsg import archive_observations
    from ESOAsg import msgs
    from ESOAsg.ancillary import cleaning_lists

    ra_list = cleaning_lists.from_element_to_list(args.ra_degree, element_type=float)
    dec_list = cleaning_lists.from_element_to_list(args.dec_degree, element_type=float)
    if args.radius is None:
        radius = None
    else:
        radius = cleaning_lists.from_element_to_list(args.radius, element_type=float)[0]
    instruments = cleaning_lists.from_element_to_list(args.instruments, element_type=str)
    data_types = cleaning_lists.from_element_to_list(args.data_types, element_type=str)

    msgs.start()
    for ra, dec in zip(ra_list, dec_list):
        position = coordinates.SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='fk5')
        msgs.newline()
        msgs.info('Query for ESO archival data around the position: {}'.format(position.to_string('hmsdms')))
        result_from_query = archive_observations.query_from_radec(position, radius=radius, instruments=instruments,
                                                                  data_types=data_types, verbose=False)
        if len(result_from_query['dp_id']) > 0:
            archive_observations.download(result_from_query['dp_id'])
    msgs.end()
