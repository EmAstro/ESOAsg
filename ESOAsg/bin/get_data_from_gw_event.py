#!/usr/bin/env python3

import argparse

from ESOAsg.ancillary import astro
from ESOAsg.core import download_archive
from ESOAsg import msgs
from ESOAsg import __version__

# from IPython import embed


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
        retrieved from the archive (default is `max_rec=9999`).
        
        A step-by-step description of the procedure is available here:
        `https://github.com/EmAstro/ESOAsg/blob/master/ESOAsg/docs/notebooks/HOWTO_getDataFromGWContours.ipynb`
                
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Input probability map out of the GW event pipeline')
    parser.add_argument('-l', '--level', nargs='+', type=float, default=50.,
                        help='Level from which extract the contours. Default is 50.')
    parser.add_argument('-s', '--show_sky', nargs='+', type=bool, default=False,
                        help='Show the contours on a sky-map')
    parser.add_argument('-a', '--asp_link', nargs='+', type=bool, default=True,
                        help='Open ASP web-pages of the selected contours')
    parser.add_argument('-d', '--download_data', nargs='+', type=bool, default=False,
                        help='Download all data collected within the considered contours')
    parser.add_argument('-i', '--instrument_name', nargs='+', type=str, default=None,
                        help='ESO instrument for which you would like to download the data')
    parser.add_argument('-m', '--max_rec', nargs='+', type=int, default=9999,
                        help='Maximum number of files you would like to download')
    parser.add_argument('-v', '--version', action='version', version=msgs._version)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        
        Create 50% confidence contours from the S191205ah_bayestar HEALPix maps, show their distribution
        on the sky, open links to the ESO ASP webpages showing data within each of the retrieved contours,
        and download the first `MUSE` cube present.
        
        get_data_from_gw_event.py S191205ah_bayestar.fits.gz -l 50. -s True -a True -d True -i MUSE -m 1
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    # File_names
    input_fits = args.input_fits
    # Confidence level
    credible_level = args.level
    # Show plot
    show_figure = args.show_sky

    msgs.start()
    contours = astro.contours_from_gw_bayestar(input_fits, credible_level=credible_level)
    astro.show_contours_from_gw_bayestar(input_fits, contours=contours,
                                        cmap='afmhot', contours_color='white', show_figure=show_figure)

    # polygons = download_archive.contours_to_polygons(contours, max_vertices=30)
    # download_archive.query_ASP_from_polygons(polygons=polygons, open_link=True)

    # results_from_TAP = download_archive.query_TAP_from_polygons(polygons=polygons,
    #                                                            merge=False, maxrec=maxrec,
    #                                                            verbose=False)
    # for idx_poly in range(0, len(results_from_TAP)):
    #     download_archive.download(results_from_TAP[0]['dp_id'])
    msgs.end()
