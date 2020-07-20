#!/usr/bin/env python3

import argparse

import numpy as np
import hashlib

from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import name_resolve
from astropy.coordinates import FK5
# from astropy import units as u

from ESOAsg import msgs
from ESOAsg.core import fitsfiles
from ESOAsg.core import download_archive

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro corrects the cubes and images of MUSE large programme 197.A-0384 to make them compatible with 
        the Phase3 standard. For more details on the project, please see
        Lofthouse et al. 2019, doi: 10.1093/mnras/stz3066
        
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
        fix_muse_p197.py --input_fits 
        """

if __name__ == '__main__':
    args = parse_arguments()

    # getting fits names
    input_fits = np.str(args.input_fits[0])
    if args.output_fits is None:
        output_fits = input_fits.replace('.fits', '_fixed.fits')
    else:
        output_fits = np.str(args.output_fits[0])

    msgs.start()

    hdul_original = fits.open(input_fits)




    msgs.work('Updating checksum and datasum')
    hdul[0].add_datasum()
    hdul[1].add_datasum()
    hdul[0].add_checksum(override_datasum=True)
    hdul[1].add_checksum(override_datasum=True)
    msgs.work('Flushing changes.')
    hdul.flush()
    hdul.close()

    msgs.newline()
    msgs.info('File {} produced.'.format(output_fits))
    msgs.end()
