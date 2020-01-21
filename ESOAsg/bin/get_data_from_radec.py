#!/usr/bin/env python3

import argparse

from astropy import coordinates
from astropy import units as u

from ESOAsg.core import download_archive
from ESOAsg import msgs

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""
        This macro collect data from the ESO archive given RA and DEC
        in degrees.

        This ueses ESOAsg version {:s}
        """.format(msgs._version),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-ra', '--ra_deg', nargs='+', type=float, default=None,
                        help='RA of the target in degree [J2000]')
    parser.add_argument('-dec', '--dec_deg', nargs='+', type=float, default=None,
                        help='Dec of the target in degee [J2000]')
    parser.add_argument('-r', '--radius', nargs='+', type=float, default=None,
                        help='Search cone radius in arcsec')
    parser.add_argument('-i', '--instrument', nargs='+', type=str, default=None,
                        help='ESO instrument')
    parser.add_argument('-v', '--version', action='version', version=msgs._version)
    return parser.parse_args()

EXAMPLES = """
        Example:
        get_data_from_radec.py --ra 15.054250 --dec 28.048833
        """

if __name__ == '__main__':
    args = parse_arguments()
    position = coordinates.SkyCoord(ra=args.ra_deg*u.degree, dec=args.dec_deg*u.degree, frame='fk5')
    
    download_archive.query_from_radec(position)
    msgs.info("End of the script")
