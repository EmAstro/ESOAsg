#!/usr/bin/env python3

import argparse
# from astropy import coordinates
# from astropy import units as u
import numpy as np
import os.path
# import glob

# from ESOAsg.core import download_archive
from ESOAsg import msgs
from ESOAsg import __version__
from ESOAsg.core import fitsfiles
# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro reads the `PROV` from a fits files and compare a selected header `card` with the ESO archive.
        
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Fits file name form which read the header. This may contain wildcards.')
    parser.add_argument('-hn', '--hdu_number', nargs='+', type=int, default=1,
                        help='select which HDU you want to check. See `fitsfiles` in `ESOAsg.core` for details.')
    parser.add_argument('-car', '--cards', nargs='+', type=str, default=None,
                        help='Header cards you want to check')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Examples:
        Check that the value of CRVAL1 in matches the one in the `PROV` files in all fits files in a folder
        check_card_from_prov.py *.fits -car CRVAL1
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    # File_names
    input_fits = args.input_fits

    # which_hdu
    from IPython import embed()
    embed()
    if isinstance(args.hdu_number, int):
        hdu_number = np.array([args.hdu_number], dtype=np.int_)
    else:
        hdu_number = np.array(args.hdu_number, dtype=np.int_)

    # cards
    if args.cards is None:
        msgs.error('At least one card should be given')
    else:
        cards = args.cards

    # Starting to read the header(s)
    msgs.start()

    for fits_file in input_fits:
        hdr = fitsfiles.header_from_fits_file(fits_file, which_hdu=hdu_number, mode='readonly')

    msgs.end()
