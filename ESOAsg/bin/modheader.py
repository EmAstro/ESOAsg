#!/usr/bin/env python3

import argparse

from astropy import coordinates
from astropy import units as u
import numpy as np
import os.path
import glob

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

    parser.add_argument('-in', '--input_fits', nargs='+', type=str, default=None,
                        help='Fits file name form which read the header. This may contain wildcards.')
    parser.add_argument('-hn', '--hdu_number', nargs='+', type=int, default=0,
                        help='select which HDU to be modified. See `fitsfiles` in `ESOAsg.core` for further details.')
    parser.add_argument('-car', '--card', nargs='+', type=str, default=None,
                        help='Header card to be modified')
    parser.add_argument('-val', '--value', nargs='+', type=str, default=None,
                        help='Value to be assigned to the card')
    parser.add_argument('-com', '--comment', nargs='+', type=str, default=None,
                        help='Comment to be assigned to the cards')
    parser.add_argument('-out', '--output', nargs='+', type=str, default=None,
                        help='Name of the output modified file')
    parser.add_argument('-ucks', '--update_checksum', nargs='+', type=bool, default=True,
                        help='Update the checksum of the file')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        mdheadercard.py -in test.fits -car CRVAL1 -val 0. -com 'This is CRVAL1' -out test_header.fits
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    if args.input_fits is None:
        msgs.error('At least one input file should be given')
    for file_name in args.input_fits:
        if not os.path.isfile(file_name):
            msgs.error('File {} does not exists'.format(file_name))
        elif not file_name.endswith('.fits'):
            msgs.error('{} needs to be a fits file'.format(file_name))
    input_fits = args.input_fits

    hdu_number = args.hdu_number

    if args.card is None:
        msgs.error('At least one card should be given')
    else:
        card = args.card

    if args.value is None:
        if len(card) == 1:
            msgs.info('The card {} will be removed'.format(card))
        else:
            msgs.info('The following cards will be removed:')
            for card_name in card:
                msgs.info(' - {}'.format(card_name))
        remove = np.bool(True)
    else:
        value = args.value
        remove = np.bool(False)
        if len(value) != len(card):
            msgs.error('There are {} values assigned to the {} cards selected'.format(len(value), len(card)))

    if args.comment is not None:
        comment = args.comment
        if len(comment) != len(card):
            msgs.error('There are {} comments assigned to the {} cards selected'.format(len(comment), len(card)))

    update_checksum = args.update_checksum

    msgs.start()


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
