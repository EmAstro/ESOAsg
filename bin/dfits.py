#!/usr/bin/env python3

import argparse

from astropy import coordinates
from astropy import units as u
import numpy as np

from ESOAsg import msgs
from ESOAsg import __version__
from ESOAsg.ancillary import cleaning_lists
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
    input_fits_files = cleaning_lists.make_list_of_fits_files(args.input_fits)

    # which_hdu
    hdu_numbers = cleaning_lists.make_list_of_int(args.hdu_number)

    # cards
    input_cards = cleaning_lists.make_list_of_strings(args.cards)

    # Starting to read the header(s)
    msgs.start()

    msgs.end()
