#!/usr/bin/env python3

import argparse
import numpy as np

from ESOAsg import msgs
from ESOAsg import __version__
from ESOAsg.ancillary import checks


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""Checking that the input fits file is valid
        
        Checks if a fits file is compliant to the FITS standard. If the option `--update_header` is enabled, the code
        will try to fix some of the header issue and the file will be modified on place.
                
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Fits file name form which read the header. This may contain wildcards.')
    parser.add_argument('-u', '--update_header', action='store_true', default=False,
                        help='Update the header on place')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        checkfits.py test.fits --update_header
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    msgs.start()

    # File_names
    for file_name in args.input_fits:
        if args.update_header:
            msgs.info('Testing and updating header: {}'.format(file_name))
            valid_fits_file = checks.fits_file_is_valid(file_name, verify_fits=True, overwrite=True)
            if not valid_fits_file:
                msgs.warning('{} not valid'.format(file_name))
                msgs.info('Header corrected')
                valid_fits_file = checks.fits_file_is_valid(file_name, verify_fits=True, overwrite=False)
        else:
            msgs.info('Testing: {}'.format(file_name))
            valid_fits_file = checks.fits_file_is_valid(file_name, verify_fits=True, overwrite=False)
        if not valid_fits_file:
            msgs.warning('{} not valid'.format(file_name))
        else:
            msgs.info('{} valid'.format(file_name))

    msgs.end()
