r"""
get_data_from_gw_event
======================
Script to collect data from the ESO archive given a `bayestar` file associated to a GW superevent

.. topic:: Inputs:

    - **input_fits** - Input probability map out of the GW event pipeline
    - **confidence_level** - Confidence level at which extract the contours
    - **max_vertices** - Max number of vertices to be considered in the conversion from contours to polygons
    - **show_sky** - Show the contours on a sky-map
    - **asp_link** - Open ASP web-pages of the selected contours
    - **download_data** - Download all data collected within the considered contours
    - **instruments** - Limit the search to the selected ESO instruments
    - **data_types** - Limit the search to the selected data type
    - **maxrec** - Maximum number of files you would like to retrieve

"""

import argparse
import matplotlib

from ESOAsg import __version__
from ESOAsg import default

STARTING_MATPLOTLIB_BACKEND = matplotlib.rcParams["backend"]


EXAMPLES = str(r"""EXAMPLES:""" + """\n""" + """\n""" +
               r"""Create the 50% confidence contours from the S191205ah_bayestar HEALPix maps, """ +
               r"""show their distribution on the sky, open links to the ESO ASP webpages showing """ +
               r"""ESO data within each of the retrieved contours, """ +
               r"""and download one of the `MUSE` cubes present: """ + """\n""" +
               r""">>> get_data_from_gw_event S191205ah_bayestar.fits.gz """ +
               r"""--confidence_level 50. """ +
               r"""--show_sky """ +
               r"""--asp_link """ +
               r"""--download_data """ +
               r"""--instruments MUSE """ +
               r"""--maxrec 1 """ + """\n""" +
               r""" """)


def parser(options=None):
    parser = argparse.ArgumentParser(
        description=r"""Collect data from the ESO archive given a `bayestar` """ +
                    r"""file associated to a GW superevent. """ + """\n""" + """\n""" +
                    r"""In summary: the code reads in a probability map in HEALPix format """ +
                    r"""(`input_fits`), finds the `confidence_level`% confidence """ +
                    r"""contours, and query the ESO archive for data within them. """ +
                    r"""Data can be showed with links to the ESO """ +
                    r"""`Archive Science Portal <http://archive.eso.org/scienceportal/home>`_ """ +
                    r"""and can be automatically downloaded on the disk. """ + """\n""" + """\n""" +
                    r"""It is possible to retrieve data for a specific instrument or data type by """ +
                    r"""setting `instruments` and/or `data_types`. """ +
                    r"""Note that `maxrec` set a limit on the number of files """ +
                    r"""that could be retrieved from the archive. """ + """\n""" + """\n""" +
                    r"""A step-by-step description of the procedure is available in this """ +
                    r"""`jupyter notebook <https://github.com/EmAstro/ESOAsg/blob/master/ESOAsg/docs/notebooks/HOWTO_getDataFromGWContours.ipynb>`_.""" + """\n""" + """\n""" +
                    r"""This uses ESOAsg version {:s}""".format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument("input_fits", nargs="+", type=str,
                        help="Input probability map out of the GW event pipeline")
    parser.add_argument("-cl", "--confidence_level", nargs="+", type=float, default=50.0,
                        help="Confidence level at which extract the contours. Default is 50.0")
    parser.add_argument("-mv", "--max_vertices", nargs="+", type=int, default=30,
                        help="Max number of vertices to be considered in the conversion from" +
                             "contours to polygons. Default it 30")
    parser.add_argument("-s", "--show_sky", action="store_true", default=False,
                        help="Show the contours on a sky-map")
    parser.add_argument("-a", "--asp_link", action="store_true", default=False,
                        help="Open ASP web-pages of the selected contours")
    parser.add_argument("-dd", "--download_data", action="store_true", default=False,
                        help="Download all data collected within the considered contours")
    parser.add_argument("-i", "--instruments", nargs="+", type=str, default=None,
                        help="Limit the search to the selected ESO instruments")
    parser.add_argument("-dt", "--data_types", nargs="+", type=str, default=None,
                        help="Limit the search to the selected data type")
    parser.add_argument("-m", "--maxrec", nargs="+", type=int, default=int(default.get_value('maxrec')),
                        help="Maximum number of files you would like to retrieve")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    return parser.parse_args()


def main(args):
    from ESOAsg.ancillary import astro
    from ESOAsg.ancillary import cleaning_lists
    from ESOAsg.ancillary import polygons
    from ESOAsg import archive_science_portal
    from ESOAsg import archive_observations
    from ESOAsg import msgs

    # Cleaning input lists
    input_fits_files = cleaning_lists.make_list_of_fits_files(args.input_fits)
    instruments = cleaning_lists.from_element_to_list(args.instruments, element_type=str)
    data_types = cleaning_lists.from_element_to_list(args.data_types, element_type=str)

    # Cleaning input values
    confidence_level = cleaning_lists.from_element_to_list(args.confidence_level, element_type=float)[0]
    maxrec = cleaning_lists.from_element_to_list(args.maxrec, element_type=int)[0]
    max_vertices = cleaning_lists.from_element_to_list(args.max_vertices, element_type=int)[0]

    # Cleaning bool
    show_figure = args.show_sky
    # Show link
    show_asp = args.asp_link
    # Download data
    download_data = args.download_data

    msgs.start()

    for fits_file in input_fits_files:
        msgs.info(' ')
        msgs.info('Working on file: {}'.format(fits_file))
        msgs.info(' ')
        contours = astro.contours_from_gw_bayestar(fits_file, credible_level=confidence_level)
        astro.show_contours_from_gw_bayestar(fits_file, contours=contours, cmap='afmhot', contours_color='white',
                                             show_figure=show_figure, matplotlib_backend=STARTING_MATPLOTLIB_BACKEND)
        polygons_list = polygons.contours_to_polygons(contours, max_vertices=max_vertices)

        if show_asp:
            msgs.info('Opening links to ASP')
            archive_science_portal.query_from_polygons(polygons=polygons_list, open_link=True)

        results_from_query = archive_observations.query_from_polygons(polygons=polygons_list, maxrec=maxrec,
                                                                      instruments=instruments, data_types=data_types,
                                                                      verbose=False)
        msgs.info(' ')
        if download_data:
            if results_from_query is None:
                msgs.warning('No data retrieved')
            elif isinstance(results_from_query, list):
                for idx_poly, result_from_query in enumerate(results_from_query):
                    msgs.info('Downloading data for polygon N.{}'.format(idx_poly+1))
                    archive_observations.download(result_from_query['dp_id'])
            else:
                msgs.info('Downloading data')
                archive_observations.download(results_from_query['dp_id'])
        msgs.info(' ')

    msgs.end()
