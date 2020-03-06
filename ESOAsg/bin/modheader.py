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
        This macro modify a card of an header and (if request) updated the value of checksum. Also, if no value is 
        given the macro just remove the card.
        
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Fits file name form which read the header. This may contain wildcards.')
    parser.add_argument('-hn', '--hdu_number', nargs='+', type=int, default=1,
                        help='select which HDU to be modified. See `fitsfiles` in `ESOAsg.core` for further details.')
    parser.add_argument('-car', '--cards', nargs='+', type=str, default=None,
                        help='Header cards to be modified')
    parser.add_argument('-val', '--values', nargs='+', type=str, default=None,
                        help='Value to be assigned to the cards')
    parser.add_argument('-com', '--comments', nargs='+', type=str, default=None,
                        help='Comment to be assigned to the cards')
    parser.add_argument('-out', '--output', nargs='+', type=str, default=None,
                        help='Name of the output modified file')
    parser.add_argument('-cks', '--checksum', nargs='+', type=bool, default=True,
                        help='Verify the checksum of the file')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Examples:
        modify the value of CRVAL1 in all the fits files in a directory and save the file putting _header.fits at the
        end.
        modheader.py *.fits -car CRVAL1 -val 0. -com 'This is CRVAL1' -out _header.fits
        
        delete the value of CRVAL1 from test0.fits and test1.fits and update the files themselves 
        modheader.py test0.fits test1.fits -car CRVAL1   
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    # File_names
    for file_name in args.input_fits:
        if not os.path.isfile(file_name):
            msgs.error('File {} does not exists'.format(file_name))
        elif not file_name.endswith('.fits'):
            msgs.error('{} needs to be a fits file'.format(file_name))
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

    # values
    if args.values is None:
        if len(cards) == 1:
            msgs.info('The card {} will be removed'.format(cards))
        else:
            msgs.info('The following cards will be removed:')
            for card in cards:
                msgs.info(' - {}'.format(card))
        remove = np.bool(True)
    else:
        values = args.values
        remove = np.bool(False)
        if len(values) != len(cards):
            msgs.error('There are {} values assigned to the {} cards selected'.format(len(values), len(cards)))

    # comments
    if args.comments is not None:
        comments = args.comments
        if len(comments) != len(cards):
            msgs.error('There are {} comments assigned to the {} cards selected'.format(len(comments), len(cards)))
    else:
        comments = [' '] * len(cards)

    # output
    if args.output is not None:
        overwrite = False
        if len(input_fits) == 1:
            output_fits = args.output
        else:
            output_fits = [output_temp.replace('.fits', np.str(args.output[0])) for output_temp in input_fits]
    else:
        output_fits = input_fits
        overwrite = True

    checksum = args.checksum

    # Starting to modify the header(s)
    msgs.start()

    for fits_file, fits_out in zip(input_fits, output_fits):
        hdul = fitsfiles.get_hdul(fits_file, mode='update', checksum=checksum)
        if remove:
            for card in cards:
                for hdu in hdu_number:
                    if card in hdul[hdu].header:
                        msgs.info('Removing header card {} from HDU N.{}'.format(card, hdu))
                        del hdul[hdu].header[card]
                    else:
                        msgs.info('Header card {} not present in HDU N.{}'.format(card, hdu))
        else:
            for card, value, comment in zip(cards, values, comments):
                for hdu in hdu_number:
                    if card in hdul[hdu].header:
                        msgs.info('Updating header card in HDU N.{}: {}'.format(hdu, card))
                        msgs.info('From {} to {} / {}'.format(hdul[hdu].header[card], value, comment))
                    else:
                        msgs.info('Adding header card in HDU N.{}: {}={} / {}'.format(hdu, card, value, comment))
                    hdul[hdu].header[card] = fitsfiles.check_value(value), comment
        # This is an astropy option to  check if the headers are broadly consistent with the standard.
        hdul.verify('fix')
        if overwrite:
            hdul.flush()
        else:
            hdul.writeto(fits_out, checksum=checksum, overwrite=True)
        hdul.close()

    msgs.end()
