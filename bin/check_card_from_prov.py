#!/usr/bin/env python3

import argparse
# from astropy import coordinates
# from astropy import units as u

from ESOAsg import msgs
from ESOAsg import __version__
from ESOAsg.core import fitsfiles
from ESOAsg.core import phase3
from ESOAsg.ancillary import cleaning_lists

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""Gather keywords from list of provenances in a fits file.""" + """\n""" + """\n""" +
                    r"""This macro reads the list of `PROVj` present in a fits file, """ +
                    r""" gets the header of the `PROVj`, and returns a selected header """ +
                    r"""`card`""" + """\n""" + """\n""" +
                    r"""This uses ESOAsg version {:s}""".format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument("input_fits", nargs="+", type=str,
                        help=r"Fits file names form which read the header. " + r"This may contain wildcards")
    parser.add_argument("-hn", "--hdu_numbers", nargs='+', type=int, default=0,
                        help=r"Select which HDU to be modified")
    parser.add_argument("-car", "--cards", nargs="+", type=str, default=None,
                        help=r"Header cards to be checked")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    return parser.parse_args()


EXAMPLES = str(r"""EXAMPLES:""" + """\n""" + """\n""" +
               r"""Check that the value of CRVAL1 in matches the one in the """ +
               r"""`PROVj` files in all fits files in a folder""" + """\n""" +
               r"""> check_card_from_prov *.fits -car CRVAL1 """ + """\n""" +
               r""" """)


if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    # File_names
    input_fits_files = cleaning_lists.make_list_of_fits_files(args.input_fits)

    # which_hdu
    hdu_numbers = cleaning_lists.make_list_of_int(args.hdu_numbers)

    # cards
    input_cards = cleaning_lists.make_list_of_strings(args.cards)

    # Starting to read the header(s)
    msgs.start()

    for fits_file in input_fits_files:
        for which_hdu in hdu_numbers:
            hdr = fitsfiles.header_from_fits_file(fits_file, which_hdu=which_hdu)
            original_files, prov_files, original_cards, prov_cards = [], [], [], []
            if 'PROV1' in hdr:
                for prov_number in hdr['PROV*']:
                    file_id = hdr[prov_number].replace('.fits', '')
                    phase3.get_header_from_archive(file_id, text_file=file_id+'.hdr')
                    hdr_prov = fitsfiles.header_from_txt_file(file_id+'.hdr')
                    for card in input_cards:
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
