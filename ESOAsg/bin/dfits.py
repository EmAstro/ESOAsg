#!/usr/bin/env python3

import argparse

from astropy import coordinates
from astropy import units as u
import numpy as np

from ESOAsg import msgs
from ESOAsg import __version__
from ESOAsg.ancillary import checks

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro mimic the behaviour of ESO's `dfits` and `grep`, i.e. to show on terminal (and on a file if
        requested) the header of a fits file. If the option `check_cards` is enabled, only the selected list of header
        cards will be printed (i.e., the result of:
        > `dfits test.fits | grep PRODCATG > test_header.txt `
        and of
        > `dfits.py test.fits -c PRODCATG -o test_header.txt`
        should be identical.
        
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Fits file name form which read the header. This may contain wildcards.')
    parser.add_argument('-hn', '--hdu_number', nargs='+', type=int, default=1,
                        help='select which HDU to be explored. See `fitsfiles` in `ESOAsg.core` for further details.')
    parser.add_argument('-car', '--cards', nargs='+', type=str, default=None,
                        help='List of header cards to be extracted')
    parser.add_argument('-o', '--output_text', nargs='+', type=str, default=None,
                        help='Result will be saved in this text file')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        dfits.py test.fits -c PRODCATG -o test_header.txt
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    # File_names
    for file_name in args.input_fits:
        if not checks.fits_file_is_valid(file_name):
            msgs.error('{} is not valid'.format(file_name))
    input_fits = args.input_fits

    # which_hdu
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

    msgs.end()
