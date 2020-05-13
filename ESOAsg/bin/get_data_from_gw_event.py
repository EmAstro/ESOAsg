#!/usr/bin/env python3

import argparse
import matplotlib

from ESOAsg.ancillary import astro
from ESOAsg.ancillary import checks
from ESOAsg.core import download_archive
from ESOAsg import msgs
from ESOAsg import __version__
from ESOAsg import default

STARTING_MATPLOTLIB_BACKEND = matplotlib.rcParams["backend"]

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro collect data from the ESO archive given a `bayestar` file associated to a GW superevent.
        
        The code reads in a probability map in HEALPix format (`input_fits`), finds the `level`% confidence
        contours, and query the ESO archive for data within them. Data can be showed with links to the ESO
        Archive Science Portal (`http://archive.eso.org/scienceportal/home`) or can be automatically 
        downloaded on the disk. 
    
        The download can be performed also for a single instrument setting `instrument_name` to a valid
        ESO instrument identifier (for instance `instrument_name=MUSE` will download all MUSE data within
        the selected contours). Note that `max_rec` set a limit on the number of files that could be 
        retrieved from the archive (default is `max_rec={}`).
        
        A step-by-step description of the procedure is available here:
        `https://github.com/EmAstro/ESOAsg/blob/master/ESOAsg/docs/notebooks/HOWTO_getDataFromGWContours.ipynb`
                
        This uses ESOAsg version {:s}
        """.format(default.get_value('maxrec'), __version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Input probability map out of the GW event pipeline')
    parser.add_argument('-l', '--level', nargs='+', type=float, default=50.,
                        help='Level from which extract the contours. Default is 50.')
    parser.add_argument('-s', '--show_sky', action='store_true', default=False,
                        help='Show the contours on a sky-map')
    parser.add_argument('-a', '--asp_link', action='store_true', default=False,
                        help='Open ASP web-pages of the selected contours')
    parser.add_argument('-d', '--download_data', action='store_true', default=False,
                        help='Download all data collected within the considered contours')
    parser.add_argument('-i', '--instrument_name', nargs='+', type=str, default=None,
                        help='ESO instrument for which you would like to download the data')
    parser.add_argument('-m', '--max_rec', nargs='+', type=int, default=default.get_value('maxrec'),
                        help='Maximum number of files you would like to download')
    parser.add_argument('-v', '--version', action='version', version=msgs._version)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        
        Create 50% confidence contours from the S191205ah_bayestar HEALPix maps, show their distribution
        on the sky, open links to the ESO ASP webpages showing data within each of the retrieved contours,
        and download the first `MUSE` cube present.
        
        get_data_from_gw_event.py S191205ah_bayestar.fits.gz -l 50. -s -a -d -i MUSE -m 1
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    # File_names
    for file_name in args.input_fits:
        if not checks.fits_file_is_valid(file_name):
            msgs.error('File {} is not a valid fits file'.format(file_name))
    input_fits_files = args.input_fits
    number_of_files = len(input_fits_files)

    # Confidence level
    if not isinstance(args.level, list):
        args.level = [args.level]
    if len(args.level) == 1 and number_of_files > 1:
        credible_levels = args.level * number_of_files
    elif len(args.level) == number_of_files:
        credible_levels = args.level
    else:
        msgs.error('A `level` should be defined for each file in input')

    # max rec
    if not isinstance(args.max_rec, list):
        args.max_rec = [args.max_rec]
    if len(args.max_rec) == 1 and number_of_files > 1:
        max_records = args.max_rec * number_of_files
    elif len(args.max_rec) == number_of_files:
        max_records = args.max_rec * number_of_files
    else:
        max_records = [default.get_value('maxrec')] * number_of_files

    # instrument
    if not isinstance(args.instrument_name, list):
        args.instrument_name = [args.instrument_name]
    if len(args.instrument_name) == 1 and number_of_files > 1:
        instruments = args.instrument_name * number_of_files
    elif len(args.instrument_name) == number_of_files:
        instruments = args.instrument_name
    else:
        msgs.error('An `instrument` should be defined for each file in input')

    # Show plot
    show_figure = args.show_sky
    # Show link
    show_asp = args.asp_link
    # Download data
    start_download = args.download_data

    msgs.start()

    for fits_file, credible_level, maxrec, instrument in zip(input_fits_files, credible_levels, max_records,
                                                             instruments):
        msgs.info(' ')
        msgs.info('Working on file: {}'.format(fits_file))
        contours = astro.contours_from_gw_bayestar(fits_file, credible_level=credible_level)
        astro.show_contours_from_gw_bayestar(fits_file, contours=contours,
                                            cmap='afmhot', contours_color='white', show_figure=show_figure,
                                            matplotlib_backend=STARTING_MATPLOTLIB_BACKEND)
        polygons = download_archive.contours_to_polygons(contours, max_vertices=30)

        if show_asp:
            msgs.info('Opening links to ASP')
            download_archive.query_ASP_from_polygons(polygons=polygons, open_link=True)

        results_from_TAP = download_archive.query_TAP_from_polygons(polygons=polygons,
                                                                    merge=False, maxrec=maxrec,
                                                                    verbose=False, instrument=instrument)
        if start_download:
            for idx_poly in range(0, len(results_from_TAP)):
                msgs.info('Downloading data for polygon N.{}'.format(idx_poly))
                download_archive.download(results_from_TAP[idx_poly]['dp_id'])

    msgs.end()
