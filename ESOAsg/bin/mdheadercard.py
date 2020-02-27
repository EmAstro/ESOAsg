#!/usr/bin/env python3

import argparse

from astropy import coordinates
from astropy import units as u
# import numpy as np

from ESOAsg.core import download_archive
from ESOAsg import msgs
from ESOAsg import __version__

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro modify a card of an header and (if request) updated the value of checksum.
        
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-i', '--input_fits', nargs='+', type=str, default=None,
                        help='Fits file name form which read the header')
    parser.add_argument('-hn', '--hdu_number', nargs='+', type=int, default=0,
                        help='select which HDU to be modified. See `fitsfiles` in `ESOAsg.core` for further details.')
    parser.add_argument('-c', '--cards', nargs='+', type=str, default=None,
                        help='List of header cards to be modified')
    parser.add_argument('-val', '--values', nargs='+', type=str, default=None,
                        help='List of values to be assigned to the cards')
    parser.add_argument('-com', '--comments', nargs='+', type=str, default=None,
                        help='List of the comments to be assigned to the cards')
    parser.add_argument('-o', '--output_text', nargs='+', type=str, default=None,
                        help='Result will be saved in this text file')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        mdheadercard.py -i test.fits -c ['CRVAL1'] -val ['0.'] -com ['This is CRVAL1'] -o test_header.txt
        """

if __name__ == '__main__':
    msgs.start()
    # input arguments and tests
    args = parse_arguments()
    input_fits = np.str(args.input_fits[0])
    if input_fits is None:
        msgs.error('At least one input file should be given')
    msgs.end()




'''    msgs.info('RA and Dec query for ESO archival data')
    msgs.newline()
    position = coordinates.SkyCoord(ra=args.ra_deg*u.degree, dec=args.dec_deg*u.degree, frame='fk5')
    result_from_query = download_archive.query_from_radec(position, args.radius)
    if args.instrument_name is not None:
        if len(args.instrument_name) > 1:
            msgs.error('Too many instrument. Only one allowed')
        instrument_name = str(args.instrument_name[0])
        msgs.info('Limit search to {} only data'.format(instrument_name))
        select_by_instrument = (result_from_query['instrument_name'].data == instrument_name.encode('ascii'))
        download_archive.download(result_from_query['dp_id'][select_by_instrument])
    if len(result_from_query['dp_id']) > 0:
        download_archive.download(result_from_query['dp_id'])
'''
