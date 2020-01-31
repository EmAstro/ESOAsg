#!/usr/bin/env python3

import argparse

import numpy as np

from astropy.io import fits
from astropy import coordinates
from astropy import units as u

from ESOAsg import msgs
from ESOAsg.core import fitsfiles
from ESOAsg.core import download_archive

from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro make the spectra of standard stars collected with MUSE compatible with the Phase3 standard.
        For more details on the project, please see:
        `https://www.eso.org/sci/publications/messenger/archive/no.178-dec19/messenger-no178-17-18.pdf`_      

        This uses ESOAsg version {:s}
        """.format(msgs._version),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-i', '--input_fits', nargs='+', type=str, default=None,
                        help='Original fits file')
    parser.add_argument('-o', '--output_fits', nargs='+', type=str, default=None,
                        help='Output fits file with corrected header and structure')
    parser.add_argument('-v', '--version', action='version', version=msgs._version)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        fix_muse_std.py --input_fits 
        """

if __name__ == '__main__':
    args = parse_arguments()

    # getting fits names
    input_fits = np.str(args.input_fits[0])
    if args.output_fits is None:
        output_fits = input_fits.replace(".fits", "_fixed.fits")
    else:
        output_fits = np.str(args.output_fits[0])

    msgs.start()

    # Copy relevant information from input into output file
    fitsfiles.new_fits_like(input_fits, [0], output_fits, overwrite=True)

    # This file will be modified in place
    hdul = fits.open(output_fits, 'update', checksum=True)
    hdr0 = hdul[0].header
    hdr1 = hdul[1].header

    msgs.work('Flushing changes.')
    hdul.flush()
    hdul.close()

    msgs.newline()
    msgs.info('File {} produced.'.format(output_fits))
    msgs.end()
