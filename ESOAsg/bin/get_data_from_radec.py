#!/usr/bin/env python3

import argparse

from astropy import coordinates
from astropy import units as u
# import numpy as np

from ESOAsg.core import archive_observations
from ESOAsg import msgs

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro collect data from the ESO archive given RA and Dec in degrees.
        
        ..note::
            ToDo EMA
            This is still work in progress, the main issue is that, when providing a negative value for Dec, `argparse` 
            has some troubles. The (ugly) workaround is to pass the argument as --dec_deg=-15.00
        
        This uses ESOAsg version {:s}
        """.format(msgs._version),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-ra', '--ra_deg', nargs='+', type=float, default=None,
                        help='RA of the target in degree [J2000]')
    parser.add_argument('-dec', '--dec_deg', nargs='+', type=float, default=None,
                        help='Dec of the target in degree [J2000]')
    parser.add_argument('-r', '--radius', nargs='+', type=float, default=None,
                        help='Search cone radius in arcsec')
    parser.add_argument('-i', '--instrument_name', nargs='+', type=str, default=None,
                        help='ESO instrument')
    parser.add_argument('-v', '--version', action='version', version=msgs._version)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        get_data_from_radec.py --ra_deg 15.054250 --dec_deg 28.048833
        get_data_from_radec.py --ra_deg="15.054250" --dec_deg=-28.048833
        """

if __name__ == '__main__':
    args = parse_arguments()
    # Check on radius
    if args.radius is None:
        radius = None
    else:
        radius = args.radius[0]
    msgs.start()
    msgs.info('RA and Dec query for ESO archival data')
    msgs.newline()
    position = coordinates.SkyCoord(ra=args.ra_deg*u.degree, dec=args.dec_deg*u.degree, frame='fk5')
    result_from_query = archive_observations.query_from_radec(position, radius)[0]
    if args.instrument_name is not None:
        if len(args.instrument_name) > 1:
            msgs.error('Too many instrument. Only one allowed')
        instrument_name = str(args.instrument_name[0])
        msgs.info('Limit search to {} only data'.format(instrument_name))
        select_by_instrument = (result_from_query['instrument_name'].data == instrument_name.encode('ascii'))
        archive_observations.download(result_from_query['dp_id'][select_by_instrument])
    if len(result_from_query['dp_id']) > 0:
        archive_observations.download(result_from_query['dp_id'])
    msgs.end()
