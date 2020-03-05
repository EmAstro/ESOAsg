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
from ESOAsg.core import download_archive

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro reads the `PROV` from a fits files and compare a selected header `card` with the ESO archive.
        
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Fits file name form which read the header. This may contain wildcards.')
    parser.add_argument('-hn', '--hdu_number', nargs='+', type=int, default=0,
                        help='select which HDU you want to check. See `fitsfiles` in `ESOAsg.core` for details.')
    parser.add_argument('-car', '--cards', nargs='+', type=str, default=None,
                        help='Header cards you want to check')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Examples:
        Check that the value of CRVAL1 in matches the one in the `PROV` files in all fits files in a folder
        check_card_from_prov.py *.fits -car CRVAL1
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    # File_names
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

    for fits_file in input_fits:
        for which_hdu in hdu_number:
            hdr = fitsfiles.header_from_fits_file(fits_file, which_hdu=which_hdu, mode='denywrite')
            original_files, prov_files, original_cards, prov_cards = [], [], [], []
            if 'PROV1' in hdr:
                for prov_number in hdr['PROV*']:
                    file_id = hdr[prov_number].replace('.fits', '')
                    download_archive.get_header_from_archive(file_id, text_file=file_id+'.hdr')
                    hdr_prov = fitsfiles.header_from_txt_file(file_id+'.hdr')
                    for card in cards:
                        original_files.append(fits_file)
                        prov_files.append(file_id)
                        original_cards.append(card)
                        prov_cards.append( hdr_prov[card])
    msgs.info('Summary:')
    msgs.info('Input <- Provenance   -  card -> Value')
    for original_file, prov_file, original_card, prov_card in zip(original_files, prov_files, original_cards,
                                                                      prov_cards):
        print('{} <- {}  -  {} -> {}'.format(original_file, prov_file, original_card, prov_card))

    msgs.end()
